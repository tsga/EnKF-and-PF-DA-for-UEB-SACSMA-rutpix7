//**********************************************************************************************
//
//  Copyright (C) 2012  David Tarboton, Utah State University, dtarb@usu.edu.  http://www.engineering.usu.edu/dtarb
//
//  This file is part of UEB.
//
//    UEB is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    UEB is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    A copy of the GNU General Public License is included in the file gpl.txt.
//    This is also available at: http://www.gnu.org/licenses/.
//
//  ifyou wish to use or incorporate this program (or parts of it) into 
//  other software that does not meet the GNU General Public License 
//  conditions contact the author to request permission.
//  David G. Tarboton  
//  Utah State University 
//  8200 Old Main Hill 
//  Logan, UT 84322-8200 
//  USA 
//  http://www.engineering.usu.edu/dtarb/ 
//  email:  dtarb@usu.edu 
//
//********************************************************************************************** 
//
#ifndef UEBDASRPFDAFUNCTIONS_H
#define UEBDASRPFDAFUNCTIONS_H

#include <ctime>
#include "Eigen/Dense"
#include<iomanip>
#include<iostream>
#include <fstream>
#include <cstring>
#include<sstream>
//#include <stdio.h>
//#include "mpi.h"
#include <netcdf.h>
//#include <netcdf_par.h>
#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
#include <random>

using namespace Eigen;

#include "eigenmvn.h"
//#include <cuda_runtime.h>
//#include "device_launch_parameters.h"
#pragma warning(disable : 4996)
//#include <algorithm>
//using namespace std;
//error handling for netcdf
#define  ERR(e) {std::cout<<"Error: "<< nc_strerror(e)<<std::endl; return 2; }
//error handling for cuda
//define  cuda_checkERR(err) { if (err != cudaSuccess) std::cout << "Error: "<<cudaGetErrorString(err)<< std::endl; exit(EXIT_FAILURE); }
//__device__ __host__
//void cuda_checkERR(cudaError_t err);

struct dasitevar {
	float elevation;
	float slope;
	float aspect;
	float cc;
	float hcan;
	float lai;
	float lat;
	float lon;
};

struct ncdaOutput {
	char outfName[256];
	char symbol[256];
	char units[256];
	int outvarIndex;
};

class uebEnKFDA { 

	public:	
		//default contr.
		uebEnKFDA();
		uebEnKFDA(int modGridCells, std::vector<std::pair<int, int> > icellCoordinates, float iy0, float ix0, const char* daconFile,
			const char* indaforcFile, const char* indaQFile, const char* uebDAsacrutpix7State);
		//uebEnKFDA(const uebEnKFDA& uebEnKFDACell0);
		//uebEnKFDA& operator= (const uebEnKFDA& uebEnKFDACell0);
		~uebEnKFDA() 
		{
			//clean up
		};

		// check whether to assimilate at current time step
		//bool daAssimlate;
		//bool updateDaArray;

		//da setttings
		//int mod_gridSize, 
		float y0, x0;
		double daTime;
		int ns_statSize;       // for now 1 state assimilated
		int numObsPoints;
		int mod_gridSize;     //model domain number of grids
		int tot_gridSize;	  // total number of grid cells includeing obs points
		int es_enseSize;	
		int es_pfSize_Fact;         // particle size = enssize * es_pfSize_Fact
		float dyC;
		float dxC;
		float forcEnStdev;    //forcing ensemble standard deviation except temperature
		float tempEnStdev;    //temperature forcing ensemble standard deviation 
		float forcCorLength;     //correlation length for forcing
		float dastateCorLength;     //states correlation length
		float tdecorrLength;     //temporal decorrelation length (hrs, dt units)
		float obsErrStdev;		  //observed (state or equivalent) var Standard deviation	
		float obsQErrStdev;		  //for Q: observed (state or equivalent) var Standard deviation	

		float daStatesStdev;      //initial model state stddev
		float daStatesStdev2;   //stdev for second (sac-rx7) state
		int daderivedType;
		int nRecs;
		//for ens run with historical forcing
		int ModelStartDateEns[4];
		int ModelEndDateEns[4];
		int forecastDateTimeEns[4];
		double forecastDateTime;   //forecast time--usually April 1
		char xmrg1dEnsDir[256];    // directory for 1d xmrg historical forcing to be used in ensemble--can be same as out dir
		int qUpdateFreq;        // Q update frequency for lagged Q assimn
		bool obsOutsideWS;		// Observations outside watershed domain used

		float dasnIndx;
		float polyThreshold;
		float polyCoeff1[5];
		float polyCoeff2[5];
		char pointInputFile[256];
	
		//int startIndexDA; 
		//int startIndexDAQ; 
		//int ncReadStartDA; // , tEndDA;		
		
		int stateIndex;    //the state to assimilate: based on array below
		char* uebDAsacrutpix7outStates[9] = { "Us", "SWE", "tausn", "Wc",  "refDepth", "totalRefDepth", "Tave", "TSURFs", "SWIT" };// 

		std::vector<ncdaOutput> daOutArr;
		std::vector<ncdaOutput> daEnsArr;
		//std::vector<std::vector<float> > daRegArray; 
		std::vector<Eigen::VectorXf > daRegArray;

		Eigen::VectorXf  daQstrArray;
		std::vector<double> daTcorrArrQ;

		std::vector<float> daYcorrArr;
		std::vector<float> daXcorrArr;
		std::vector<double> daTcorrArr;

		//debug file
		//std::ofstream debugOutputFile;
		Eigen::VectorXi Hc_hgVector;  //Hc= grid cell indices with observation
		void setHcMatrices(std::vector<std::pair<int, int> > icellCoordinates);
		
		void readDaContr(const char* daconFile);
		void readDAOutputControl(const char* outputconFile);
		void  readTStextFile_multiVal(const char* inforcFile);   // , int &numdaPoints);    // , daControlV &inpDaControlV);
		void readTStextFileQ(const char* inforcFile);
		/*__host__ __device__*/
		void updateDaArr(int& startIndexDA);
		/*__host__ __device__*/
		void initDAMatrices(std::vector<std::pair<int, int> > cellCoordinates);		
		void initDAMatrices_Default(std::vector<std::pair<int, int> > cellCoordinates);
		/*__host__ __device__*/
		void runEnKF(Eigen::VectorXf Z_obs, std::vector<Eigen::RowVectorXf> Xh_obsState, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> ensObservationErr, Eigen::RowVectorXf stateOutputArr, Eigen::RowVectorXf& stateOutputUpdate);       // , bool NormalDist);  //float* &ensembleUpdateArr,
		//11.23.18 each obs points da without its observation (Leave one approach)
		/*__host__ __device__*/
		void runEnKF_LeaveOneOut(int igo, Eigen::VectorXf Z_obs, std::vector<Eigen::RowVectorXf> Xh_obsState, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> ensObservationErr, Eigen::RowVectorXf stateOutputArr, Eigen::RowVectorXf& stateOutputUpdate);       // , bool NormalDist);  //float* &ensembleUpdateArr,
		//9.6.18: for obs points
		void runEnKFopt(int ido, Eigen::VectorXf Z_obs, Eigen::RowVectorXf Xh_obsState, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> ensObservationErr, Eigen::RowVectorXf stateOutputArr, Eigen::RowVectorXf& stateOutputUpdate);         // bool NormalDist)   //float* &ensembleUpdateArr,
		
		int julian(int yy, int mm, int dd);
		double julian(int I, int M, int K, double H);
		// flag to write warnings,...etc 
		//int uebDAsacrutpix7_outflag;    //output print
		static const int uebDAsacrutpix7_debugout1 = 0;  // print major debug information
		static const int uebDAsacrutpix7_debugout2 = 0;  // more detailed print 
		static const int uebDAsacrutpix7_debugout3 = 0;  // print in 'main' 
		
		//private:
		Eigen::EigenMultivariateNormal<float> std_norm_dist_Forc_Default;
		Eigen::EigenMultivariateNormal<float> norm_dist_0Mean_Tempr;     //for temperature
		Eigen::EigenMultivariateNormal<float> norm_dist_0Mean_Default;
		Eigen::EigenMultivariateNormal<float> norm_dist_1Mean_Default;
		Eigen::EigenMultivariateNormal<float> norm_dist_SACRX;  // for sac states

		Eigen::VectorXf Z_obs;
		Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> R_obsErrCov;
		//Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> P_stateCov;
		float P_stateCov;
		//std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > P_stateCovBackground;
		//std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > P_stateCovBackground_Points;   //10.11.17 at obs points
};
#endif


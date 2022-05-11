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
#ifndef UEBDASRPFUEBPGDECLS_H
#define UEBDASRPFUEBPGDECLS_H

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
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#pragma warning(disable : 4996)

using namespace Eigen;

//error handling for netcdf
#define  ERR(e) {std::cout<<"Error: "<< nc_strerror(e)<<std::endl; return 2; }
#define  ERR2(e, FILE_NAME) {std::cout<<"Error: "<< nc_strerror(e)<< FILE_NAME<< std::endl; std::getchar(); return 2; }
//error handling for cuda
//define  cuda_checkERR(err) { if (err != cudaSuccess) std::cout << "Error: "<<cudaGetErrorString(err)<< std::endl; exit(EXIT_FAILURE); }

//41618 for sacsma in gpu

/*__host__ __device__
void fland1Device(float *pxv, float *edmnd, float *ta, float *we, float *aesc, float *sh, float *dtday,
	float *sacst, float *frzst, float *sacpar, float *frzpar, int *nsoil, int *nupl, int *nsac, int *ivers,
	float *surf, float *grnd, float *tet, float *smc, float *sh2o, float *sacst_prv__,
	float *dtfrz, int *idtfrz, float *frzdup, float *frzdbt, float *frost, float *tsint, float *swint, float *swhint,
	float *dsint, int *ndsint, float *dsintw, int *ndintw,
	float *sif, float *bfp, float *bfs, float *ssur, float *sdro, float *sperc,
	int *normalize, float *hrapx, float *hrapy, int *error);*/

__host__ __device__
void fland1Device(float pxv, float dt, float aesc,
	float sacst[6], float sacpar[17],
	float &edmnd, float &surf, float &grnd, float &tet,
	float &sif, float &bfp, float &bfs, float &ssur, float &sdro, float &sperc,
	float hrapx, float hrapy, int &error);

class SACGridDA {
	public:
		SACGridDA() {
			//const
		}
		~SACGridDA() {
			//dest
		}
		float pxv;
		float dtm, dtday;
		float aesc;

		float sacst[6] = { 0 };
		float sacpar[17];
		float edmnd = 0;
		float surf = 0;
		float grnd = 0;
		float tet = 0;
		//float edmndOrg = edmnd;

		float sif, bfp, bfs, ssur, sdro, sperc;
		float hrapx, hrapy;
		int error;

		// for ens simulation
		float sacstEns[300][6];   // , sacst2[300], sacst3[300], sacst4[300], sacst5[300], sacst6[300];
		float surfEns[300] = { 0 };
		float grndEns[300] = { 0 };

		__host__ __device__
		void runSACGrid_Ens(int grdIndx, int nEns, float *pxvEns);

		__host__ __device__
		void runSACGrid()
		{
			//dtday = dtm / (24. * 60);  /* convert time step in days */
			//aesc = 0;
			fland1Device(pxv, dtday, aesc,
				sacst, sacpar, 
				edmnd, surf, grnd, tet,
				sif, bfp, bfs, ssur, sdro, sperc,
				hrapx, hrapy, error);
		}
		__host__ __device__
		void updateSACGrid(std::vector< float > &sac_st, float &totalWater1,	float &totalWater2);
		
};

int do_sac_EnsHost_Ens(int threadsPerBlock, int npix, int nEns, SACGridDA *sacGridArray, float *pxvEns);

int do_sac_EnsHost(int threadsPerBlock, int npix, SACGridDA *sacGridArray);

void do_sac_EnsHost_PC(int npix, SACGridDA *sacGridArray);


struct params {
	float irad, ireadalb, tr, ts, ems, cg, z, zo, rho, rhog, lc, ks, de, avo,
	anir0, lans, lang, wlf, rd1, dnews, emc, alpha, alphal, gpar, uc, as, Bs,
	lambda, rimax, wcoeff, apar, cpar;
};

struct sitevar {
	char svName[256];
	int   svType;
	char svFile[256];
	char svVarName[256];
	float svdefValue;
	float** svArrayValues;
};

struct inpforcvar {
	char infName[256];
	int  infType;
	char infFile[256];
	char infvarName[256];
	char inftimeVar[256];
	float infdefValue;
	int numNcfiles;
};

struct pointOutput {
	char outfName[256];
	int ycoord;
	int xcoord;
};

struct ncOutput {
	char outfName[256];
	char symbol[256];
	char units[256];
};

struct aggOutput {
	//char outfName[256];
	char symbol[256];
	char units[256];
	char aggop[256];
};

struct inptimeseries {
	//CTime dtime;
	float datetime;
	float tsValue;
};

class uebCellDA {  

	public:
		uebCellDA();
		uebCellDA(float Params[32], int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT,
			double  UTCoffset, int inpDailyorSubd, int oStride);
		//uebCellDA(uebCellDA& uCell0);
		//uebCellDA& operator= (uebCellDA& uCell0);
		~uebCellDA();

		//for EnKF
		float forcEnStdev;    //forcing ensemble standard deviation except temperature
		float tempEnStdev;
		float corrFact1;   // = 1.0 - modelDT / decorLength;
		float corrFact2;	// = sqrtf(1.0 - corrFact1 * corrFact1);    //1.0 - corrFact1;  // 
		float tdecorLength;
		//
		int startIndexDA;
		//for PF
		int startIndexDAQ;
		// check whether to assimilate at current time step
		bool daAssimlate;
		bool updateDaArray;
		bool updateDaArrayQ;
		//bool unitFahrenheit;
		float modisAlbedoFact;  // albedo multiplier for modis

	   //read parameter values
		__host__ __device__
		void  setConstantValues();
		void  readParams(const char* inpFile, float Params[32]);
		__host__ __device__
		void  setParams(float Params[32]);
		//copy site variable array
		__host__ __device__
		void  setSiteVars_and_Initconds(float SiteVars[32]);
		__host__ __device__
		void  setModelRun_Settings(int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT, 
			                       double  UTCoffset, int inpDailyorSubd, int outtStride);
		__host__ __device__
		void setInitialEnsembleStates(int nEns);		//, const char* forcName);
		  //##### snowdv.cpp
		__host__ __device__
		int getforcOffset(int ioffst, int dimLen2);
		 __host__ __device__
		 void setForcingAndRadiationParamterization();	
         //11.27.17 for RDHM only one time step considered at a time
		 //inputs are P, Ta, Tmin(daily), Tmax(daily), VP (humidity), V
		 __host__ __device__
		 void setForcingAndRadiationParamterization(float P, float Ta, float Tmin, float Tmax, float VP, float V);	
		 /*__host__ __device__*/ 
		 //save only state with observation
		 //void runUEBEnsembles(int nEns, int outStateIndex, float* &stateOutput);
		 __host__ __device__
		 void runUEB();
		 //when historical forcing--including rad parametrization
		 __host__ __device__
		 void runUEBEnsembles_HistForcing(int nEns);
		 // 3.20.18 for ens run on device side: copy spatially correlated random samples
		 __host__ __device__
		void copyEnsForcingMultiplier(int nEns, std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > ensForcingMultiplier, int cellRank);
		 // 8.27.18 --corr among forc--for ens run on device side : copy spatially correlated random samples
		__host__ __device__
		void copyEnsForcingMultiplier(int cellRank, int nEns, int totalGridLength, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>  ensForcingMultiplier,
			std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > ensForcingMultiplierVRH);
		__host__ __device__
		void copyEnsForcingMultiplier(int cellRank, int nEns, int totalGridLength, float *dev_multivarNormalDistSamplesForc,
				float *dev_multivarNormalDistSamplesForcTVRH0, float *dev_multivarNormalDistSamplesForcTVRH1, float *dev_multivarNormalDistSamplesForcTVRH2);
		 // 3.20.18 for ens run on device side: copy historical forc
		 __host__ __device__
		 void copyEnsForcingHist(float *uebForcEnsArray, int ensNum, int cellRank, int numCells);

		 //3.20.18 ens random multipliers copied as class members
		 __host__ __device__
		void runUEBEnsembles(int nEns);
		 /*__host__ __device__*/ 
		 //1.10.18 for only one state (eg. SWE)
		 //void updateBackgroundStates(int nEns, float*** ensForcingMultiplier, int outStateIndex, int cellRank, float* updateStateArr);
		 //set background states for next time step to filter updated
		 __host__ __device__
		 void updateBackgroundStates(int nEns, int outStateIndex, Eigen::RowVectorXf updateStateArr);
		 //4.11.18 update by weight
		 __host__ __device__
		void updateBackgroundStates(int nEns, int outStateIndex, Eigen::RowVectorXf updateStateArr, std::vector<int> weightIndices);
		 
		 __host__ __device__
		 void setNextStepStates(int nEns);
		 /*__host__ __device__*/ 
		 //void runUEB(int dimLen2);	//this uses class forcing arrays	 
		 __host__ __device__ 
		 void updateSimTime();

		 //11.27.17--for RDHM single time step considered at a time 
		 /*__host__ __device__*/
		 //void updateSimTime(int Year, int Month, int Day, double dHour);

		 void readInputForContr(const char* inputconFile);		
		 //void getInpForcArr(int numNc[13], float*** RegArray[13], float &tcorVar, int ncTotaltimestep[13], MPI::Intracomm inpComm, MPI::Info inpInfo);
		 void updateInpForcArr(float*** RegArray[13], int ncTotaltimestep[13]);
		 void setInpForcArr(int it, float ***inArray, float* forcArr, int ncTotaltimestepit);
         // read input forcing time series text file
		 void  readTStextFileTimeValPair(const char* inforcFile, std::vector<std::pair<double, float> > &tvar_in, int &nrecords);
         //### 
		 void  printPointOutputs(const char* outFileName);
		 void  printDebugOutputs();
		 void  printSampleOutputs(const char* outFileName);

		//accumulation zone?
		bool accumulationZone;
		int modelStartDate[3], modelEndDate[3];
		double modelStartHour, modelEndHour, modelDT, UTCOffset, modelSpan;
		//1.18.18 moved to public
		double currentModelDateTime;

		//int ModelStartDate[3], ModelEndDate[3], TstepsinDay, numTotalTimeStepss; //check this 
		//double ModelStartHour, ModelEndHour, ModeldT, ModelUTCOffset, modelSpan;
		//track current grid cell
		int uebCellDAX;
		int uebCellDAY;
		//Global variables ###_TBC 5.6.13
		// flag to write warnings,...etc ###***TBC 9.20.13
		int snowdgtvariteflag;  // 0; print debug information
		int snowdgtvariteflag2;  // more detailed print in Surfebsc
		int snowdgtvariteflag3;  // print in predictor corrector function
		int snowdgt_outflag;        //output print
		int radwarnflag;
		int outtStride;
		int numTotalTimeSteps;
		int numSimTimeSteps;           // number of steps to simulated before updating data from netcdf (24 except towards the end of simulation)
		int istep;                   //simulation step, 0...23 for a batch of simualations inbetween data update
		int nstepinaDay;
		//int inpDailyorSubdaily;       // 0: values given at each (sub-daily time steps); 1: daily values
		//int numSimTimeSteps;
		int timeSeriesIndex;       //indicates whether a time series is read or not --- forcing time series from text file need to be read only once
		//float* tsvarArray[13]; 
		inpforcvar infrContArr[13]; 		
		//int ncTotaltimestep[13];
		int startIndex[13], ncReadStart[13], tEnd;		
		//outputs
		float OutVarValues[71];                 //11.30.17 71st Trange    //4.28.15 compute max of 96 time steps at time (for daily input, 4 days with dT = 6h)
		//params and site and inp vars
		float statesiteValues[32], paramValues[32];

		//historical forcing for ensemble run
		static const int ensSize = 60;
		float  PrecArr[ensSize], TempArr[ensSize], TaminArr[ensSize], TamaxArr[ensSize], WindspArr[ensSize], RhArr[ensSize],
			VpArr[ensSize], ApresArr[ensSize], SradArr[ensSize], LradArr[ensSize], NradArr[ensSize], QgArr[ensSize], SnowalbArr[ensSize];
		
		float  tsprevday[24], taveprevday[24];   //assumes the dt >= 1hr which is generally the case in UEB
		//time related funcs
		//snowxv.cpp
		//********UPDATEtime ()  Update time for each time step
		 __host__ __device__  
		void  UPDATEtime(int &YEAR, int &MONTH, int &DAY, double &HOUR, double DT);
		//    function to return number of days in February checking for leap years
		__host__ __device__ 
		int lyear(int year);		
		//To convert the real date to julian date
		__host__ __device__  
		int julian(int yy, int mm, int dd);		
		//these were copied from functions.f90
		//COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  
		__host__ __device__ 
	    double julian(int I, int M, int K, double H);
		//COMPUTES CALENDAR DATE AND time, GIVEN JULIAN DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
		__host__ __device__  
		void  calendardate(double TJD, int &I, int &M, int &K, double &H);

		//std::vector<float> stateVS[9];    ////9.17.17 9th for RMelt=SWIT   // 0, stateV1, stateV2, stateV3, stateV4, stateV5;    //8.8.16 7th state for snow surface temp, 8th Tave
		float stateVS[9][60];   //3.20.18 for EnKF on device side
		float ensForcingMultiplier[6][60];

    private:
		//to keep states not updated
		float stateVS0[9];   ////9.17.17 9th for RMelt=SWIT 
		//std::vector<float> cumPV, cumEsV, cumEcV, cumMrV, cumGmV, cumEgV, dStorageV, errMBV;
		float cumPV[60], cumEsV[60], cumEcV[60], cumMrV[60], cumGmV[60], cumEgV[60], dStorageV[60], errMBV[60];
		//int inputforcIndex;    // the forcing to perturb for ensembles--using the string lists below
		char* uebForcInpts[8] = { "Ta", "P", "V", "RH", "Qsi", "Qli", "Qnet", "cosZen" };   //"atff", "HRI", "Eacl", "Ema"

	    //variables to track previous day temperature profile
		//std::vector <float> tsprevday;
		//std::vector <float> taveprevday;
		//js  Constant    floatset 
		float T_0; 				// Temperature of freezing (0 C)
		float T_k;  // 273.15;			// Temperature to convert C to K (273.15)
		float SB_c;  // Stefan boltzman constant (2.041334e-7 KJ/m^2-hr-K^4) #corrected 12.23.14
		float H_f;  // 333.5;			// Heat of fusion (333.5 KJ;  // kg)
		float Hne_u;  // 2834.0;		// Heat of Vaporization (Ice to Vapor, 2834 KJ;  // kg)
		float C_w;  // 4.18;			// Water Heat Capacity (4.18 KJ;  // kg;  // C)
		float C_s;  // 2.09;			// Ice heat capacity (2.09 KJ;  // kg;  // C)
		float C_p;  // 1.005;			// Air Heat Capacity (1.005 KJ;  // kg;  // K)
		float Ra_g;  // 287.0;			// Ideal Gas constant for dry air (287 J;  // kg;  // K)
		float K_vc;  // 0.4;			// Von Karmans constant (0.4)
		float Hs_f;  // 3600.0;		// Factor to convert ;  // s into ;  // hr (3600)
		float Rho_i;  // 917.0;			// Density of Ice (917 kg;  // m^3)
		float Rho_w;  // 1000.0;		// Density of Water (1000 kg;  // m^3)
		float Gra_v;  // 9.81;			// Gravitational acceleration (9.81 m;  // s^2)
		float W1da_y;  // 0.261799;		// Daily frequency (2pi;  // 24 hr 0.261799 radians;  // hr) 
		float  Io;  // 4914.0;            //  Solar constant  Kj/m^2/hr
		//pi copied from snowdxv.f90
		float P_i;  // 3.141592653589793238462643383279502884197169399375105820974944592308;		// Pi
		//data for pred-corr
		float wtol;  // 0.025;
		float utol;  // 2000.0;
		//from TURBFLUX()
		float tol;  // 0.001;
		int  nitermax;  // 20;  
		int  ncitermax; //21      
		//common (surface and average temp of previous time step ? ### 9.20.13?)
		float Tsk_save;  // 273.16,
		float Tssk_old;  // 273.16,  
		float Tsavek_old;  // 273.16, 
		float Tsavek_ave;  // 273.16, 
		float Tssk_ave;  // 273.16/*added 6.7.13*/;
//==============================================changes for new conf 5.1.15
//from snowdv
		//const int nsv = 10, npar = 32, nxv = 6, niv = 8, nov = 14;
		float  Inpt[8], sitev[10], outv[14], statev[6], Param[32], dtbar[12], mtime[4];                                     // YJS pass to reflect the change of Snow (Year, month, Date, Hour)
		int irad, ireadalb, subtype, iradfl, iflag[6];
		float  slope, azi, bca, bcc, lat, ts_last, lon;
		//float *tsprevday, *taveprevday;
		//int stepinaDay, nstepinaDay; 
		//float referenceHour, referenceTime, CTJD;
		double dHour, EJD, MHour, sHour, UTCHour, OHour, UTCJulDat;  //currentModelDateTime, 1.18.18 moved to public
		float fHour, fModeldt;
		int Year, Month, Day, MYear, MMonth, MDay;

		// CHANGES TO ACCOMODATE GLACIER
		float  WGT; // WGT=WATER EQUIVALENT GLACIER THICKNESS
		float Us, Ws, Wc, Apr, cg, rhog, de, tave, Ws1, Wc1, cumP, cumEs, cumEc, cumMr, cumGm, cumEg, Ta, P, V, RH, Tmin, Tmax, Trange, Qsiobs, Qg, Qli, QLif;
		float Vp; //#12.18.14 Air vapor pressure 
		float Qnetob, Snowalb, cosZen, atff, cf, atfimplied, Ema, Eacl, dStorage, errMB, HRI0, HRI, as, bs;
		float  OutArr[53];
		//to get max/min daily temperatures
		int nb, nf, nbstart, nfend;
		double dstHour;	

		//snowxv.cpp		
		//    to get the atmospheric transmissivity using the Bristow and Campbell  (1984) approach
		__host__ __device__  
		void  atf(float &atff,float trange,int month, float *dtbar, float a, float c);
		// To get hourly radiation index
		__host__ __device__  
	    void  hyri(int YEAR, int MONTH, int DAY, float HOUR, float DT, float SLOPE, float AZI, float LAT, float &HRI, float &COSZEN);		 
		//    Computes the incoming longwave radiation using satterlund Formula
		__host__ __device__  
		void  cloud(float as, float bs, float atff, float &cf);
		//???? long wave radiation from temperatrue and other weather variables??
		//TBC_6.5.13
		__host__ __device__  
		void  qlif(float TA, float RH, float TK, float SBC, float &Ema, float &Eacl, float cf, float &qliff );	

		//##### canopy.cpp       
		//     Partition the incoming solar radiation into direct and diffuse components
		__host__ __device__  
	    void  PSOLRAD( float Qsi, float atff, float *param, float cf,     
		                float &Taufb, float &Taufd, float &Qsib, float &Qsid);     // Output variables
		//     Estimates the direct and diffuse solar radiation fractions transmitted through the canopy      
		__host__ __device__  
		void  TRANSRADCAN (float COSZEN, float *sitev, float *param,       
		     	            float &Betab, float &Betad, float &Taub, float &Taud);                     // Output variables: Transmission and Reflection fractions
		//      Computes the net canopy and sub-canopy solar radiation
		//      considering the multiple scatterings of radiation by the canopy and multile reflections 
		//      between the canopy and the snow surface
		__host__ __device__  
		void  NETSOLRAD(float Ta,float A,float Betab,float Betad,float Wc,float Taub,float Taud, 
					     float Qsib,float Qsid,float *param, float	Fs,float &Qsns,float &Qsnc); //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
		//     Computes the net canopy and beneath canopy longwave radiation
		__host__ __device__ 
		void  NETLONGRAD(float RH,float Ta,float Tss,float Tc,float Tk,float Fs,float EmC,float EmS,float SBC,float cf,float *sitev,float Qli,float *param,
						    float &Qlis, float &Qlns,float &Qlnc );
		//****************  BULK AERODYNAMIC AND CANOPY BOUNDARY LAYER RESISTANCES  *******************
		//      Calculates resistances for the above and beneath canopy turbulent heat fluxes 
		__host__ __device__  
		void  AeroRes(float P, float Wc, float V, float Ta, float Tss, float Tc, float Fs, float *param, float *sitev, float Tk,
					    float &d, float &Z0c, float &Vz, float &Rc, float &Ra, float &Rbc, float &Rl, float &RKINc, float &RKINa, float &RKINbc, float &RKINl);         // Output variables
		//     Calculates wind at the 2 m above the surface and at the sink using the input of measured wind at 2 m above the canopy
		__host__ __device__  
		void  WINDTRANab(float V, float Zm, float *param, float *sitev,float &Vstar,float &Vh,float &d,float &Z0c,float &Vz,float &Vc);               // Out put wind speed within canopy at height 2 m
		//     Calculates the turbulent heat fluxes (sensible and latent heat fluxes) and condensation/sublimation.
		__host__ __device__  
		void  TURBFLUX (float Ws, float Wc, float A, float Tk,float Tc,float Ta,float Tss, float RH, float V,float Ea,float P,float *param,float *sitev,
		                  float &d, float &Z0c, float &Vz, float &Rkinc, float &Rkinbc, float &Tac, float &Fs, float &Ess, float &Esc,  // Output variables				 
		                  float &QHc, float &QEc, float &Ec, float &QHs, float &QEs, float &Es, float &QH, float &QE, float &E ); 
		//     Calculates amount of snow intercepted by the canopy
		__host__ __device__  
		void  INTERCEPT (float Ta, float LAI,float P, float Wc, float dt, float Inmax, float Uc, float Cc,
		     	         // Output variables
						 float &ieff, float &Ur, float &intc
					   );
		//     Calculates the heat advected to the snowpack due to rain
		__host__ __device__  
		float QPF(float PR, float TA, float TO, float PS, float RHOW, float HF, float CW, float CS);
		//      Routine to correct energy and mass fluxes when numerical overshoots dictate that W was changed in the calling routine - either because W went negative
		//      or due to the liquid fraction being held constant.
		__host__ __device__  
		void  PREHELP( float W1,float W,float DT,float &FM,float FM1,float fac,float PS,float PRAIN,float &E,float RHOW,float HF,float &Q,float &QM,float &MR, float &QE, float HSF);
		//     Computes the exponential integral function for the given value      
		__host__ __device__ 
		float EXPINT(float LAI);
		//     Calculates the vapour pressure at a specified temperature over water or ice depending upon temperature.  Temperature is celsius here.
		__host__ __device__
		float svp(float T);
		//     Calculates the vapour pressure at a specified temperature over water using polynomial from Lowe (1977).
		__host__ __device__ 
		float svpw(float T);
		//     Calculates the vapour pressure at a specified temperature over ice using polynomial from Lowe (1977).
		__host__ __device__ 
		float svpi(float T);
		//     Estimates reflection and scattering coefficient
		__host__ __device__ 
		float Tau1(float Rho, float G, float h, float COSZEN, float kk);
		__host__ __device__ 
		float Tau2(float Rho, float G, float h, float kk, float EXPI);

		//#######  snowdgtv.cpp
		//ueb 'driving' function
		//yjs Note: in this subroutine, the outv is an array which passes value to this subroutine and back to the snow 
		//yjs drive program. The Outv(9) and outv(10) pass the daily average snowpack temperature and daily
		//yjs snow surface temperature to this subroutine but pass the Qin total and combined mass fluxes 
		//yjs back. //
		__host__ __device__  
		void  SNOWUEB2(); 
		    //float dt, float *input, float *sitev, float *statev, float* tsprevday, float* taveprevday, int &nstepday, float *param, int *iflag,
			//float &cump, float &cumes, float &cumEc, float &cumMr, float &cumGM, float &cumEg, float *outv, float *mtime, float atff, float cf, float *OutArr);
		//************************  Predictor CORRECTOR ***********************************
		__host__ __device__  
		void  PREDICORRc ( float &Us, float &Ws, float &Wc, float Alb, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float Qli,
						    float atff, float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime, 
		          	        // Following variables are output
		                   float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qmc, float &Mc, float &FMc,float &intc, float &Inmax, float &ieff, float &Ur, float &Cf, 
						float &Taufb, float &Taufd, float &Qsib, float &Qsid, float &Taub,float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &Tac,
						// Just for testing
					       float &QHs, float &QEs,float &Es, float &QPs, float &MR, float &QMs, float &Q,float &FM, float &TSURFs, float &tave, float &Qnet, float &refDepth, float &totalRefDepth, float &smelt,
						   float &gsurf, float &Qlis
					   );
		//				CALCULATE CANOPY MASS FLUX AT ANY INSTANT
	    __host__ __device__  
		void  QFMc(float Us, float Ws, float Wc, float A, float dt, float rid, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float Qli,
				      float cosZen, float EmC, float Ems, float *param, float *sitev, int iradfl, float Qnetob, int iTsMethod, float *mtime,
					   // Following variables are output
					     float &Tc, float &QHc, float &QEc, float &Ec, float &Qpc, float &Qps, float &Qmc, float &Mc, float &FMc, float &intc,float &Inmax, float &ieff, float &Ur, float &Cf, float &Taufb,
						 float &Taufd, float &Qsib, float &Qsid, float &Taub, float &Taud, float &Qsnc, float &Qsns, float &Qlnc, float &Qlns, float &Rkinc, float &Rkinsc, float &Vz, float &TSURFs, float &Tac,
						 // Just for testing
						 float &tave, float &qnet, float &QHs, float &QEs, float &Es, float &MR, float &QMs, float &Q, float &FM, float &refDepth, float &totalRefDepth, float &Smelt, float &smeltc
				 );
		//  COMPUTE THE SURFACE TEMPERATURE OF SNOW 
		__host__ __device__ 
		float SRFTMPSC(float &Tssk, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float cf, float Qli, float cosZen,
						   float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Qpsin,float Inmax,float Rkinc,float Rkinsc,float Vz,
						     float &Tc,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &Tave,float &refDepth,float &smelt,float &smeltC
					    );
		// COMPUTE THE CANOPY TEMPERATURE
		__host__ __device__ 
		float CanTemp(float& Tck, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Qsi, float atff, float cf, float Qli,
		               float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpcin,float Inmax,float Rkinsc,float Vz,
					     float& Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &Fs,float &tave,float &refDepth,float &smeltC
					);
		//      FUNCTION TO EVALUATE THE SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR SURFACE TEMPERATURE                      DGT and C Luce 4/23/97
		__host__ __device__ 
		float SURFEBSC(float Tssk, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Fs, float Cf, float Qli, float Qsi,
		                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float Qpc,float Qps,float &Inmax,
						    float &Rkinc,float &Rkinsc,float &Vz,float Tck,float Tk,float &Tak,float &EA,float &RHOA,float &fkappaS,float &RHO,float &TherC,float &TSURFs,float &Tave,float &refDepth
						);// Heat and vapor conductance for neutral
		//*************  FUNCTION TO EVALUATE THE CANPPY SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR CANOPY TEMPERATURE 
		__host__ __device__ 
		float SURFEBC(float Tck, float Us, float Ws, float Wc, float A, float dt, float P, float Pr, float Ps, float Ta, float V, float RH, float Fs, float Cf, float Qli, float Qsi,
		                  float atff,float cosZen,float EmC,float Ems,float *param,float *sitev,int iradfl,float Qnetob,int iTsMethod,float *mtime,float QPc,float &QPs,float &Inmax,
						     float &Rkinc,float &Rkinsc,float &Vz,float Tssk,float Tk,float &Tak,float &EA,float &RHOA,float &FkappaS,float RHO, float &TherC,float &TSURFs,float &tave,float &refDepth
					  );                             // Reduced parametre later                                    // Heat and vapor conductance for neutral
		//                           ********  QcEst() *******
		__host__ __device__
		float QcEst(float Ws, float P, float Tssk, float &Tck, float V, float &Zm, float &d, float &Z0c, float &Rimax, float &Rcastar, float Cf, float Fs, float Qli,
		               float &Hcan, float Vz,float Ta, float Rh, float RKINsc, float QPs, float To, float Ps,float Qsi, float atff, float cosZen, float APr,float Tak,
					     float &EA, float &A, float &Ac, float Wc, float Inmax, float Qnetob, int Iradfl, float *param, float *sitev
					);
		//     Calculates the melt rate and melt outflow
		__host__ __device__ 
		float FMELT(float UB, float RHOW, float W, float HF, float LC, float RID, float KS, float PRAIN);
		// Linear simple function to calculate the gradient of Q and T
		__host__ __device__  
		void  Grad(float qc1, float qc2,float t1,float t2, float &a,float &b);
		//     function to calculate the refreezing depth
		__host__ __device__ 
		float  refDep(float flans, float a, float b, float hf, float  rhom, float dt, float x1);
		__host__ __device__ 
		float TAVG(float UB, float W, float RHOW, float CS, float TO, float RHOG, float DE, float CG, float HF);
		//yjs  Calculate the daily average value //yjs  n number of records, a minimum value -100 or some
		__host__ __device__ 
		float daily_ave(float* backup, int n, float a);
		//	Function to get the LanE which is the thermal conductivity by ze
		__host__ __device__ 
		float  LanE(float LanS, float LanG, float Zs, float rho, float rhog, float cs, float cg, float r, float &ze, float w1day);
		//     Function to calculate Albedo  BATS Albedo Model (Dickinson et. al P.21)
		__host__ __device__ 
		float ALBEDO(float tausn, float coszen, float d, float aep, float abg, float avo, float airo);
		//     Function to calculate Dimensionless age of snow for use in BATS Albedo Model (Dickinson et. al P.21)
		__host__ __device__  
		void  AGESN(float &tausn, float dt, float Ps, float tsurf, float tk, float dNewS);
		//     Partitioning of precipitation into rain and snow      
		__host__ __device__ 
		float PARTSNOW(float P, float TA, float TR, float TS);



		//utility functions
		//findMin findMax functions
		__host__ __device__
			float findMax(float a, float b);
		__host__ __device__
			float findMin(float a, float b);
		__host__ __device__
			int findMax(int a, int b);
		__host__ __device__
			int findMin(int a, int b);
		//print array values (for inspection during debugging)
		/*// __host__ __device__  void  printArrayvalues(int arrLength, float* arrValues);
		// __host__ __device__  void  printArrayvalues(int yLength, int xLength, float **arrValues);*/
};
//perturbed forcing used upto forecast date
int callUEBrunHostEnsemble_PertForcing(int threadsPerBlock, int npix, int nEns, uebCellDA *uebCellDAObjArr, float *uebForcArray);
//for hist forcing ens
int callUEBrunHostEnsemble_HistForcing(int threadsPerBlock, int npix, int nEns, uebCellDA *uebCellDAObjArr, int nOpts, uebCellDA *daUebCellObjArr);  // , float* uebForcArrayDaPoint);  //,uebEnKFDA objEnKFArr          
//3.25.18 calll gupfrom host (obs points separately simulated)
int callUEBrunHostEnsemble(int threadsPerBlock, int npix, int nOpts, int nEns, uebCellDA *uebCellDAObjArr, uebCellDA *daUebCellObjArr, float *uebForcArray, float* uebForcArrayDaPoint, 
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>  ensForcingMultiplier, std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > ensForcingMultiplierVRH);  //,uebEnKFDA objEnKFArr    
// whne no da
int callUEBrunHost(int threadsPerBlock, int npix, uebCellDA *uebCellDAObjArr, float *uebForcArray);  //,uebEnKFDA objEnKFArr 
//2.25.18 writes vector to a 2D nc at time step
__host__ __device__
int Write1DVector_to2DNC(const char* FileName, const char* varName, int tIndx, int Nz_dim, float* var_inp); //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//--without-mpi //1.25.18 writes a 2D array at time step
__host__ __device__
int Write2DSlub_to3DNC(const char* FileName, const char* varName, int tIndx, int Np_dim, int Nz_dim, float** var_inp); //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//no-mpi: 2.25.18 for ens and da outputs  //creates 2D netcdf and stores dimension variables; called once for a given output netcdf 
__host__ __device__
int Create2DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
	const char* zName, int tDim, int zDim, float* t_inp, float* fillVal);           // MPI::Intracomm inpComm, MPI::Info inpInfo)
//no-mpi: 1.18.18 for ens and da outputs //creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
__host__ __device__
int Create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
	const char* yxName, const char* zName, int tDim, int yxDim, int zDim, float* t_inp, float* fillVal);
// 10.6.17 for reading point site variable arrays from nc
__host__ __device__
int read_Point_SiteVars_NC(const char* FILE_NAME, const char* varName, float* &pvar_in); //, MPI::Intracomm inpComm, MPI::Info inpInfo)
// 10.6.17 for reading point site variable arrays---all at once-- from nc
__host__ __device__
int read_Point_SiteVars_NC(const char* FILE_NAME, float** &pvar_in); //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//1.18.18 for forcing at obs points
//read wole matrix at once
__host__ __device__
int read2DNC_Contigious(const char* FILE_NAME, const char* VAR_NAME, float** &pvar_in);  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//10.6.17 read multiple arrays (vectors of vals at points) for given time range
__host__ __device__
int readNC_vector(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, int tEnd, float** &pvar_in, int &nrecords);
//2.19.18 read array at an index (time step)
__host__ __device__
int readNC_Array_atIndex(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int tIndex, float* pvar_in);  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
//for single point //2.19.18 read array at an index (time step)
__host__ __device__
int readNC_Array_atIndices(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int tIndex, int pIndex, float* &pvar_in);
/*
//read multiple slubs (y,x arrays) along the time dim (for interval between tstart --- tend) 
// __host__ __device__ 
int readNC_yxSlub(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, int tEnd,float*** &pvar_in, int &nrecords, int &numNc, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read multiple slubs (y,x arrays) along the time dim for data with t, y, x config---time as slowely varying array; and pvar_in already allocated
// __host__ __device__ 
int readNC_yxSlub_givenT(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, float** &pvar_in, float &tcorvar, int &numNc, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read whole nc array for a given variable to contiguous 1d array; the y,x,time order is irrelevant inside the function, but it is assumed known by caller,
// __host__ __device__ 
int readNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &Ntdim, float* &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read whole nc array for a given variable to contiguous 3d array; the y,x,time order is irrelevant inside the function, but it is assumed known by caller,
// __host__ __device__ 
int read3DNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &Ntdim, float*** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read nc to contiguous 1d array;  y,x,time coordinate names are same as index names
// __host__ __device__ 
int readNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* tcorvar, float* ycorvar, float* xcorvar, size_t &Ntdim, size_t &Nydim, size_t &Nxdim, float* &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//read 3d nc to contiguous array;  y,x,time coordinate names are same as index names
// __host__ __device__ 
int read3DNC_Contiguous(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* tcorvar, float* ycorvar, float* xcorvar, size_t &Ntdim, size_t &Nydim, size_t &Nxdim, float*** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//writes the 1D array (TS values) at specified location in the netcdf
// __host__ __device__
int WriteTSto3DNC(const char* FileName, const char* VarName, int dimOrd, int y_dim, int x_dim, int Nt_dim, float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//writes multiple 2D arrays (e.g., TS-Ens values) at specified locations in the netcdf
// __host__ __device__ 
int WriteTStoMultiDnc_Block(const char* FileName, const char* VarName, int dimOrd, int *YindArr, int *XindArr, int bSize, int Nt_dim, int Nz_dim, float*** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);  
//writes multiple 1D arrays (TS values) at specified locations in the netcdf
// __host__ __device__ 
int WriteTSto3DNC_Block(const char* FileName, const char* VarName, int dimOrd, int *YindArr, int *XindArr, int bSize, int Nt_dim, float** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//aggregate outputs at specified zone in the netcdf
// __host__ __device__ 
int Write_uebaggTS_toNC(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp);
//aggregate outputs at specified zone in the netcdf in parallel 
// __host__ __device__ 
int Write_uebaggTS_toNC_par(const char* FileName, const char* VarName, int dimOrd, int z_dim, int Nt_dim, float* var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 3D netcdf in PARALLEL mode, and stores dimension variables for UEB aggregated outputs; some attributes are copied from the watershed netCDF
// __host__ __device__ 
int create3DNC_uebAggregatedOutputs(const char* FileName, aggOutput *aggOut, int naggOut, const char* tName, const char* tUnits, const char* tlong_name, const char* tcalendar, int Nt_dim, int dimOrd,
	float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, int nZones, const char * zName, float* y_inp, float* x_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 3D netcdf and stores dimension variables for a given UEB output; called once for a given output netcdf 
//attributes are copied from the 2D (watershed) netCDF which the 3D netCDF spatial dimension follow
// __host__ __device__ 
int create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnits, const char* tlong_name,
	const char* tcalendar, int Nt_dim, int dimOrd, float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName, const char* yName, const char* xName, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates a multiD netcdf to store data assimilation outputs for UEB; called once for a given output netcdf 
// __host__ __device__ 
int createMultiDnc_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnits,
	const char* tlong_name, const char* tcalendar, int Nt_dim, int dimOrd, float* t_inp, float *fillVal, const char* ws_FileName, const char* ws_VarName,
	const char* yName, const char* xName, const char* zName, int Nz_dim, MPI::Intracomm inpComm, MPI::Info inpInfo);
//creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
// __host__ __device__ 
int Create3DNC(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* yName, const char* xName, const char* tUnits,
	const char* yUnits, const char* xUnits, int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float* fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo);
//This one creates and stores array passed to it at once
// __host__ __device__ 
int Write3DNC(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* yName, const char* xName,
	const char* tUnits, const char* yUnits, const char* xUnits,
	int Nt_dim, int Ny_dim, int Nx_dim, int dimOrd, float* t_inp, float* y_inp, float* x_inp, float*** var_inp, MPI::Intracomm inpComm, MPI::Info inpInfo);
// __host__ __device__ 
int read3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME,
	const char* ycor_NAME, const char* tcor_NAME, float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//function to read multiple blocks of single column/rod along time dimension from 3D netcdf file, for given y , x coordinate arrays
// __host__ __device__ 
int readNC_TS_Block(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float** &pvar_in, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo, int *YindArr, int *XindArr, int bSize);
//function to read single column/rod along time dimension from 3D netcdf file, for given y , x coordinates
// __host__ __device__ 
int readNC_TS(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, const char* ycor_NAME, const char* xcor_NAME,
	float* &pvar_in, float* &tcorvar, int ydim, int xdim, int &nrecords, MPI::Intracomm inpComm, MPI::Info inpInfo);
// __host__ __device__ 
int read2DNC(const char* FILE_NAME, const char* VAR_NAME, float** &pvar_in, MPI::Intracomm inpComm, MPI::Info inpInfo);
//the following is to read watershed file
// __host__ __device__ 
int readwsncFile(const char* FILE_NAME, const char* VAR_NAME, const char* ycor_NAME,
	const char* xcor_NAME, float* &ycorvar, float* &xcorvar, int** &pvar_in, int &ydim, int &xdim, int &fillVal, MPI::Intracomm inpComm, MPI::Info inpInfo);
// __host__ __device__ 
int Write3DNC(const char* FILE_NAME, const char* VAR_NAME, const char* xcor_NAME,
	const char* ycor_NAME, const char* tcor_NAME, float* xcorvar, float* ycorvar, float* tcorvar, float*** pvar_out, MPI::Intracomm inpComm, MPI::Info inpInfo);

*/
//ueb inputs 
 void  readParams(const char* inpFile, float* &parArray, const int nParams);
 void  readParams(const char* inpFile, params strParamValues);
 void  readSiteVars(const char* inpFile, sitevar svArr[], int indx);
//another overload because the code fails at run time for the above_#7.4.13
 void  readSiteVars(const char* inpFile, sitevar *&svArr);
 void  readSiteVars(const char* inpFile, float svSValue[], char* svFile[], char* svVarName[], int svType[] );
 void  readInputForcVars(const char* inputconFile, inpforcvar *frArray);
 void readOutputControl(const char* outputconFile, pointOutput* &pOut, int &npout, ncOutput* &ncOut, int &nncout, aggOutput* &aggOut, int &naggOut, ncOutput* &daOut, int &ndaout);
// // __host__ __device__  void  readTextData(const char* inforcFile, inptimeseries* &strinpts, int &nrecords);
 void  readTextData(const char* inforcFile, float *&tcor_var, float *&tvar_in, int &nrecords);  
 void  readTextData(const char* inforcFile, float *&tvar_in, int &nrecords);
 void  readTStextFile(const char* inforcFile, float *&tvar_in, int &nrecords);

//arrays and matrices
__host__ __device__ 
float*** Create3DArray(int nt, int nr, int nc);
__host__ __device__  
void  Delete3DArray(float ***A, int nt, int nr, int nc);
//create 2D array and allocate contiguous memory block this enbales a full block read of netcdf
__host__ __device__ 
float** create2DArray_Contiguous(int nr, int nc);
//delets a 2D array (frees memory allocated) contiguously
__host__ __device__  
void  delete2DArray_Contiguous(float** myMatrix);
//create 3D array and allocate contiguous memory block this enbales a full block read of netcdf
__host__ __device__ 
float*** create3DArrayblock_Contiguous(int nt, int nr, int nc);      //inputs: no. of rows and no. of colos , height/time (dimensions) of matrix
//delets a 3D array (frees memory) allocated contiguously
__host__ __device__  
void  delete3DArrayblock_Contiguous(float*** myMatrix); // int nr, int nc)// input: 3D array

#endif
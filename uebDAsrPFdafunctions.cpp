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

#include "uebDAsrPFdafunctions.h"
#include <ctime>

uebEnKFDA::uebEnKFDA(int modGridCells, std::vector<std::pair<int, int> > icellCoordinates, float iy0, float ix0, const char* daconFile,
	const char* indaforcFile, const char* indaQFile, const char* uebDAsacrutpix7State)
{
	//daAssimlate = true;
	//updateDaArray = true;

	y0 = iy0;
	x0 = ix0;
	daTime = 99999.0;
	//1.13.18 TBCL: this is 1 for now 
	ns_statSize = 1;   //=inpDaControl.ns_statSize;
					   //1.12.18 augument the state vector with point obs states
					   //mod_gridSize = cellCoordinates.size();
	readDaContr(daconFile);

	forecastDateTime = julian(forecastDateTimeEns[0], forecastDateTimeEns[1], forecastDateTimeEns[2], (double)forecastDateTimeEns[3]);

	mod_gridSize = modGridCells;
	tot_gridSize = mod_gridSize + numObsPoints;    // mod_gridSize + numObsPoints;

	//startIndexDA = 0;
	//startIndexDAQ = 0;
	//ncReadStartDA = 0;
	nRecs = 1;    // inpDaControl.nRecs;
				  //8.8.16 :  //9.15.17 SWIT not state but added for ESP
	stateIndex = 1;
	for (int vindx = 0; vindx < 9; vindx++)
	{
		if (strcmp(uebDAsacrutpix7State, uebDAsacrutpix7outStates[vindx]) == 0) { //TODO: try to do without looping?
			stateIndex = vindx;
			break;
		}
	}
	//only call initDAMatrices after reading obs file
	readTStextFile_multiVal(indaforcFile);
	//for Q assmn
	readTStextFileQ(indaQFile);

	initDAMatrices(icellCoordinates);			//, daYcorrArr, daXcorrArr, uebDAsacrutpix7State);
	
	std::ofstream debugOutputFile;
	debugOutputFile.open("debugOutput.txt", std::ios::out);
	debugOutputFile.close();
}

uebEnKFDA::uebEnKFDA()
{
	//daAssimlate = true;
	//updateDaArray = true;
	y0 = 0.0;
	x0 = 0.0;
	daTime = 99999.0;
	//1.13.18 TBCL: this is 1 for now 
	ns_statSize = 1;

	ModelStartDateEns[0] = 2008;
	ModelStartDateEns[1] = 10;
	ModelStartDateEns[2] = 1;
	ModelStartDateEns[3] = 7;
	ModelEndDateEns[0] = 2009;
	ModelEndDateEns[1] = 10;
	ModelEndDateEns[2] = 1;
	ModelEndDateEns[3] = 7;

	forecastDateTime = julian(ModelStartDateEns[0], ModelStartDateEns[1], ModelStartDateEns[2], (double)ModelStartDateEns[3]);   //forecast time--usually April 1
	std::strcpy(xmrg1dEnsDir, "./outdaHF");
	obsOutsideWS = false;
	//1.12.18 augument the state vector with point obs states
	//mod_gridSize = cellCoordinates.size();
	mod_gridSize = 1;
	numObsPoints = 1;
	tot_gridSize = mod_gridSize + numObsPoints;
	es_enseSize = 1;
	es_pfSize_Fact = 1;
	dyC = 0.25;			//HRAP 800.0;
	dxC = 0.25;		//800.0;	
	forcEnStdev = 0.1;    //forcing ensemble standard deviation except temperature
	tempEnStdev = 1.0;    //temperature forcing ensemble standard deviation 
	forcCorLength = 16.7979;	//HRAP 80000.0;     //correlation length for forcing
	dastateCorLength = 16.7979;		// 80000.0;
	tdecorrLength = 24.0;
	obsErrStdev = 0.001;		  //observed (state or equivalent) var Standard deviation	
	obsQErrStdev = 2.0;		// cms
	qUpdateFreq = 16;

	daStatesStdev = 0.5;
	daStatesStdev2 = 0.2;
	daderivedType = 0;

	//startIndexDA = 0;
	//startIndexDAQ = 0;
	//ncReadStartDA = 0;
	nRecs = 1;
	//8.8.16 :  //9.15.17 SWIT not state but added for ESP
	stateIndex = 1;

	daTcorrArr.push_back(99999.0);
	daRegArray.resize(1);
	daRegArray[0].resize(1);
	daRegArray[0][0] = 0.0;
	daYcorrArr.push_back(99999.0);
	daXcorrArr.push_back(99999.0);
	std::vector<std::pair<int, int> > icellCoordinates;
	icellCoordinates.push_back(std::make_pair(0, 0));

	//8.28.18 --- initDAMatrices_Default(icellCoordinates);    //, daYcorrArr, daXcorrArr, uebDAsacrutpix7State);

	daTcorrArrQ.push_back(99999.0);
	daQstrArray.resize(1);
	daQstrArray(0) = 0.0;

	std::ofstream debugOutputFile;
	debugOutputFile.open("debugOutput.txt", std::ios::out);
	debugOutputFile.close();
}

//##### TODO: Consistency in coordinates !!!!!
void uebEnKFDA::initDAMatrices(std::vector<std::pair<int, int> > cellCoordinates)
{
	//set initial state covariance	
	/*Matrix<float, Dynamic, Dynamic, RowMajor> statesCorr(ns_statSize, ns_statSize);
	statesCorr.setZero();
	for (int i = 0; i < ns_statSize; i++) {
	statesCorr(i, i) = 1.0;
	}
	P_stateCov.resize(ns_statSize, ns_statSize);
	P_stateCov.setZero();
	//TODO: Revise this later
	for (int i = 0; i < ns_statSize; i++) {
	for (int j = 0; j < ns_statSize; j++) {
	P_stateCov(i, j) = statesCorr(i, j) * daStatesStdev * daStatesStdev;  // daStatesStdev[i] * daStatesStdev[j];
	}
	}*/
	P_stateCov = daStatesStdev * daStatesStdev;
	//std::cout << std::endl << "intitial covariance matrix" << std::endl;
	//std::cout << P_stateCov << std::endl;

	/*for (int irc = 0; irc < npix; irc++)
	P_stateCovBackground.push_back(P_stateCov);
	for (int irc = 0; irc < numObsPoints; irc++)
	P_stateCovBackground_Points.push_back(P_stateCov);    //10.11.17 at obs points
	*/

	//1.13.18 coordinates of all grid cells (obs points + sim model grids)
	std::vector<float> yCoordArray;
	std::vector<float> xCoordArray;
	yCoordArray.resize(tot_gridSize);
	xCoordArray.resize(tot_gridSize);
	int im = 0;
	for (int jn = 0; jn < cellCoordinates.size(); ++jn)
	{
		yCoordArray[im] = y0 + cellCoordinates[jn].first * dyC;
		xCoordArray[im] = x0 + cellCoordinates[jn].second * dxC;
		im++;
	}
	for (int io = 0; io < numObsPoints; io++)
	{
		yCoordArray[im] = daYcorrArr[io];
		xCoordArray[im] = daXcorrArr[io];
		im++;
	}

	//for (int iy = 0; iy < tot_gridSize; iy++) std::cout << yCoordArray[iy] << " "; std::cout << std::endl;

	// forcing covariance matrix based on distance between grid cells
	Matrix<float, Dynamic, Dynamic, RowMajor> covarFd(tot_gridSize, tot_gridSize);   // distance-based correlation
	double rij = 0;      //distance between grid cells 
	// covariance matrix based on distance between grid cells
	for (int i = 0; i < tot_gridSize; i++) {
		for (int j = 0; j < tot_gridSize; j++) {
			rij = (yCoordArray[i] - yCoordArray[j]) * (yCoordArray[i] - yCoordArray[j])			//(cellCoordinates[i].first - cellCoordinates[j].first)*(cellCoordinates[i].first - cellCoordinates[j].first) * dyC * dyC    //(i1-i2)^2 * dy^2 + (j1-j2)^2 *dx^2
				+ (xCoordArray[i] - xCoordArray[j]) * (xCoordArray[i] - xCoordArray[j]);		// (cellCoordinates[i].second - cellCoordinates[j].second)*(cellCoordinates[i].second - cellCoordinates[j].second) * dxC * dxC;
			rij = sqrt(rij);
			covarFd(i, j) = exp(-1.0 * rij / forcCorLength);							//*daContArr.forcEnStdev * daContArr.forcEnStdev;
		}
	}
	//std::cout << " Distance based correlation: " << std::endl;
	//std::cout << covarFd << std::endl;
	//std::cout << " In Dafunc 1" << std::endl;

	//for perturbation of forcing
	Matrix<float, Dynamic, Dynamic, RowMajor> covarFs(3, 3);   // correlation-between states
	covarFs <<   1.0, -0.8, 0.5,
				-0.8, 1.0, -0.5,
				 0.5, -0.5, 1.0;
/*	covarFs << 1.0000000, - 0.1018587,  0.3927145,  0.5934511,
		     - 0.1018587,  1.0000000, - 0.7979352,  0.5018560,
		       0.3927145, - 0.7979352,  1.0000000, - 0.4927249,
		       0.5934511,  0.5018560, - 0.4927249,  1.0000000;*/

	//for all forcing except temp, mean of 1 with std.dev from user input
	float forcStdDevind[6] = { tempEnStdev, forcEnStdev, forcEnStdev, forcEnStdev, forcEnStdev, forcEnStdev };
	VectorXf meanF(3 * tot_gridSize);
	meanF.setZero();
	//std::cout << " Mean vector : " << std::endl;	//std::cout << meanF << std::endl;
	Matrix<float, Dynamic, Dynamic, RowMajor> covarF(3 * tot_gridSize, 3 * tot_gridSize);
	int igrid, jgrid, istate, jstate;
	for (int i = 0; i < 3 * tot_gridSize; i++) {
		igrid = i % tot_gridSize;
		istate = i / tot_gridSize;
		for (int j = 0; j < 3 * tot_gridSize; j++) {
			jgrid = j % tot_gridSize;
			jstate = j / tot_gridSize;
			covarF(i, j) = covarFd(igrid, jgrid) * covarFs(istate, jstate);   // *forcStdDevind[istate] * forcStdDevind[jstate];
			//covar(i,j) = cov_grid(i,j) * cov_state(i,j) * Sig(i) * Sig(j);
		}
	}
	const uint64_t seedF = static_cast<uint64_t>(time(0));
	//9.1.16 set rand. generator
	//8.23.18 def seed: std_norm_dist_Forc_Default.setSeed(seedF);		
	std_norm_dist_Forc_Default.setMean(meanF);
	std_norm_dist_Forc_Default.setCovar(covarF, true);	
	//std::cout << " In Dafunc 2" << std::endl;
	
	//Ta wind and RH
	VectorXf meanF_TVRH(tot_gridSize);
	meanF_TVRH.setZero();
	//9.1.16 set rand. generator
	const uint64_t seedFVRH = static_cast<uint64_t>(time(0));	
	//8.23.18 def seed: norm_dist_1Mean_Default.setSeed(seedFVRH);		
	norm_dist_1Mean_Default.setMean(meanF_TVRH);
	norm_dist_1Mean_Default.setCovar(covarFd, true);
	//std::cout << " In Dafunc 3" << std::endl;
//
	Z_obs.resize(numObsPoints);
	R_obsErrCov.resize(numObsPoints, numObsPoints);
	//not used yet vR_obsErr.resize(mo_obseSize, es_enseSize);
	//##TODO Revise later
	R_obsErrCov.setZero(); // = covarM.cast<float>() * daContArr.obsErrStdev * daContArr.obsErrStdev;
	for (int io = 0; io < numObsPoints; io++)
	{
		R_obsErrCov(io, io) = obsErrStdev * obsErrStdev;
	}
	//for measurement / observatin 
	VectorXf meanM(numObsPoints);
	// Create a multi variate normal distribution with mean 0                    8.28.16
	meanM.setZero();  //8.28.16                             //8.28.16 for obs use y' = y + vR, vR ~ N(0,R)	
	const uint64_t seedM = static_cast<uint64_t>(time(0));
	//8.23.18 def seed: norm_dist_0Mean_Default.setSeed(seedM);
	norm_dist_0Mean_Default.setMean(meanM);
	norm_dist_0Mean_Default.setCovar(R_obsErrCov, true);
	//std::cout << " In Dafunc 4" << std::endl;

	// for sac states
	Eigen::VectorXf meanSRx7(mod_gridSize);
	meanSRx7.setOnes();
	//srx7ErrCov.setOnes();
	Matrix<float, Dynamic, Dynamic, RowMajor> covarSd(mod_gridSize, mod_gridSize);
	rij = 0;      //distance between grid cells 
    // covariance matrix based on distance between grid cells
	for (int i = 0; i < mod_gridSize; i++) {
		for (int j = 0; j < mod_gridSize; j++) {
			rij = (yCoordArray[i] - yCoordArray[j]) * (yCoordArray[i] - yCoordArray[j])			//(cellCoordinates[i].first - cellCoordinates[j].first)*(cellCoordinates[i].first - cellCoordinates[j].first) * dyC * dyC    //(i1-i2)^2 * dy^2 + (j1-j2)^2 *dx^2
				+ (xCoordArray[i] - xCoordArray[j]) * (xCoordArray[i] - xCoordArray[j]);		// (cellCoordinates[i].second - cellCoordinates[j].second)*(cellCoordinates[i].second - cellCoordinates[j].second) * dxC * dxC;
			rij = sqrt(rij);
			covarSd(i, j) = exp(-1.0 * rij /dastateCorLength);						//*daContArr.forcEnStdev * daContArr.forcEnStdev;
		}
	}
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> 
		srx7ErrCov = covarSd * daStatesStdev2 * daStatesStdev2;     //* covarSs(istate, jstate) 
	const uint64_t seedMsrx = static_cast<uint64_t>(time(0));
	//8.23.18 def seed: std_norm_dist_1D.setSeed(seedMsrx);
	norm_dist_SACRX.setMean(meanSRx7);
	norm_dist_SACRX.setCovar(srx7ErrCov, true);
	//std::cout << " In Dafunc 5" << std::endl;

	setHcMatrices(cellCoordinates);	

	return;
}
//##### TODO: Consistency in coordinates !!!!!
void uebEnKFDA::initDAMatrices_Default(std::vector<std::pair<int, int> > cellCoordinates)
{
	//set initial state covariance	
	/*Matrix<float, Dynamic, Dynamic, RowMajor> statesCorr(ns_statSize, ns_statSize);
	statesCorr.setZero();
	for (int i = 0; i < ns_statSize; i++) {
	statesCorr(i, i) = 1.0;
	}
	P_stateCov.resize(ns_statSize, ns_statSize);
	P_stateCov.setZero();
	//TODO: Revise this later
	for (int i = 0; i < ns_statSize; i++) {
	for (int j = 0; j < ns_statSize; j++) {
	P_stateCov(i, j) = statesCorr(i, j) * daStatesStdev * daStatesStdev;  // daStatesStdev[i] * daStatesStdev[j];
	}
	}*/
	P_stateCov = daStatesStdev * daStatesStdev;
	//std::cout << std::endl << "intitial covariance matrix" << std::endl;
	//std::cout << P_stateCov << std::endl;

	/*for (int irc = 0; irc < npix; irc++)
	P_stateCovBackground.push_back(P_stateCov);
	for (int irc = 0; irc < numObsPoints; irc++)
	P_stateCovBackground_Points.push_back(P_stateCov);    //10.11.17 at obs points
	*/

	//1.13.18 coordinates of all grid cells (obs points + sim model grids)
	std::vector<float> yCoordArray;
	std::vector<float> xCoordArray;
	yCoordArray.resize(tot_gridSize);
	xCoordArray.resize(tot_gridSize);
	int im = 0;
	for (int jn = 0; jn < cellCoordinates.size(); ++jn)
	{
		yCoordArray[im] = y0 + cellCoordinates[jn].first * dyC;
		xCoordArray[im] = x0 + cellCoordinates[jn].second * dxC;
		im++;
	}
	for (int io = 0; io < numObsPoints; io++)
	{
		yCoordArray[im] = daYcorrArr[io];
		xCoordArray[im] = daXcorrArr[io];
		im++;
	}

	//for (int iy = 0; iy < tot_gridSize; iy++) std::cout << yCoordArray[iy] << " "; std::cout << std::endl;

	// forcing covariance matrix based on distance between grid cells
	//for perturbation of forcing
	VectorXf meanF, meanT, meanM;
	meanF.resize(tot_gridSize);
	//for all forcing except temp, mean of 1 with std.dev from user input
	meanF.setOnes();

	Matrix<float, Dynamic, Dynamic, RowMajor> covarF;
	covarF.resize(tot_gridSize, tot_gridSize);
	double rij = 0;      //distance between grid cells 
	double corij = 1;    //correlation between grid cells i and j
						 // covariance matrix based on distance between grid cells
	for (int i = 0; i < tot_gridSize; i++) {
		for (int j = 0; j < tot_gridSize; j++) {
			rij = (yCoordArray[i] - yCoordArray[j]) * (yCoordArray[i] - yCoordArray[j])			//(cellCoordinates[i].first - cellCoordinates[j].first)*(cellCoordinates[i].first - cellCoordinates[j].first) * dyC * dyC    //(i1-i2)^2 * dy^2 + (j1-j2)^2 *dx^2
				+ (xCoordArray[i] - xCoordArray[j]) * (xCoordArray[i] - xCoordArray[j]);		// (cellCoordinates[i].second - cellCoordinates[j].second)*(cellCoordinates[i].second - cellCoordinates[j].second) * dxC * dxC;
			rij = sqrt(rij);
			corij = exp(-1.0*rij / forcCorLength);
			covarF(i, j) = corij;							//*daContArr.forcEnStdev * daContArr.forcEnStdev;
		}
	}
	Matrix<float, Dynamic, Dynamic, RowMajor> covarFS;
	covarFS = covarF * forcEnStdev * forcEnStdev;
	// Create a multi variate standard normal distribution 
	const uint64_t seedF = static_cast<uint64_t>(time(0));
	//9.1.16 set rand. generator
	//8.23.18 def seed: std_norm_dist_Forc_Default.setSeed(seedF);
	std_norm_dist_Forc_Default.setMean(meanF);
	std_norm_dist_Forc_Default.setCovar(covarFS, true);

	//for temperature use additive termwith mean 0
	meanT.resize(tot_gridSize);
	//for all forcing except temp, mean of 1 with std.dev from user input
	meanT.setZero();
	Matrix<float, Dynamic, Dynamic, RowMajor> covarTS;
	covarTS = covarF * tempEnStdev * tempEnStdev;
	// Create a multi variate standard normal distribution 
	const uint64_t seedT = static_cast<uint64_t>(time(0));
	//8.23.18 def seed: norm_dist_0Mean_Tempr.setSeed(seedT);
	norm_dist_0Mean_Tempr.setMean(meanT);
	norm_dist_0Mean_Tempr.setCovar(covarTS, true);

	Z_obs.resize(numObsPoints);
	R_obsErrCov.resize(numObsPoints, numObsPoints);
	//not used yet vR_obsErr.resize(mo_obseSize, es_enseSize);
	//##TODO Revise later
	R_obsErrCov.setZero(); // = covarM.cast<float>() * daContArr.obsErrStdev * daContArr.obsErrStdev;
	for (int io = 0; io < numObsPoints; io++)
	{
		R_obsErrCov(io, io) = obsErrStdev * obsErrStdev;
	}
	//for measurement / observatin 
	meanM.resize(numObsPoints);
	// Create a multi variate normal distribution with mean 0                    8.28.16
	meanM.setZero();  //8.28.16                             //8.28.16 for obs use y' = y + vR, vR ~ N(0,R)	
	const uint64_t seedM = static_cast<uint64_t>(time(0));
	//8.23.18 def seed: norm_dist_0Mean_Default.setSeed(seedM);
	norm_dist_0Mean_Default.setMean(meanM);
	norm_dist_0Mean_Default.setCovar(R_obsErrCov, true);

	//for perturbing sac-rx7 states Eigen::EigenMultivariateNormal<double> std_norm_dist_1D;  // for sac states
	Eigen::VectorXf meanSRx7;
	meanSRx7.resize(mod_gridSize);
	meanSRx7.setOnes();
	//srx7ErrCov.setOnes();
	Matrix<float, Dynamic, Dynamic, RowMajor> covarSt;
	covarSt.resize(mod_gridSize, mod_gridSize);
	rij = 0;      //distance between grid cells 
	corij = 1;    //correlation between grid cells i and j
				  // covariance matrix based on distance between grid cells
	for (int i = 0; i < mod_gridSize; i++) {
		for (int j = 0; j < mod_gridSize; j++) {
			rij = (yCoordArray[i] - yCoordArray[j]) * (yCoordArray[i] - yCoordArray[j])			//(cellCoordinates[i].first - cellCoordinates[j].first)*(cellCoordinates[i].first - cellCoordinates[j].first) * dyC * dyC    //(i1-i2)^2 * dy^2 + (j1-j2)^2 *dx^2
				+ (xCoordArray[i] - xCoordArray[j]) * (xCoordArray[i] - xCoordArray[j]);		// (cellCoordinates[i].second - cellCoordinates[j].second)*(cellCoordinates[i].second - cellCoordinates[j].second) * dxC * dxC;
			rij = sqrt(rij);
			corij = exp(-1.0*rij / dastateCorLength);
			covarSt(i, j) = corij;							//*daContArr.forcEnStdev * daContArr.forcEnStdev;
		}
	}
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> srx7ErrCov;
	srx7ErrCov = covarSt * daStatesStdev2 * daStatesStdev2;
	const uint64_t seedMsrx = static_cast<uint64_t>(time(0));
	//8.23.18 def seed: std_norm_dist_1D.setSeed(seedMsrx);
	norm_dist_SACRX.setMean(meanSRx7);
	norm_dist_SACRX.setCovar(srx7ErrCov, true);


	setHcMatrices(cellCoordinates);

	return;
}

//this finds the indices of the grid cells where there are observations
void uebEnKFDA::setHcMatrices(std::vector<std::pair<int, int> > icellCoordinates)
{
	//if (useHmatrix) 
	Hc_hgVector.resize(numObsPoints);
	Hc_hgVector.setConstant(mod_gridSize / 2);    //middle of the array as default ()
	int im = 0;
	//point observations
	double rij = 0.0;
	double minDist1 = (icellCoordinates[mod_gridSize - 1].first * dyC) * (icellCoordinates[mod_gridSize - 1].first * dyC)			//(cellCoordinates[i].first - cellCoordinates[j].first)*(cellCoordinates[i].first - cellCoordinates[j].first) * dyC * dyC    //(i1-i2)^2 * dy^2 + (j1-j2)^2 *dx^2
		+ (icellCoordinates[mod_gridSize - 1].second * dxC) * (icellCoordinates[mod_gridSize - 1].second * dxC);		// (cellCoordinates[i].second - cellCoordinates[j].second)*(cellCoordinates[i].second - cellCoordinates[j].second) * dxC * dxC;
	minDist1 = sqrt(minDist1);
	for (int ida = 0; ida < numObsPoints; ++ida)
	{
		double minDist = minDist1;
		for (int jn = 0; jn < mod_gridSize; ++jn)
		{
			rij = (daYcorrArr[ida] - (y0 + icellCoordinates[jn].first * dyC)) * (daYcorrArr[ida] - (y0 + icellCoordinates[jn].first * dyC))			//(cellCoordinates[i].first - cellCoordinates[j].first)*(cellCoordinates[i].first - cellCoordinates[j].first) * dyC * dyC    //(i1-i2)^2 * dy^2 + (j1-j2)^2 *dx^2
				+ (daXcorrArr[ida] - (x0 + icellCoordinates[jn].second * dxC)) * (daXcorrArr[ida] - (x0 + icellCoordinates[jn].second * dxC));		// (cellCoordinates[i].second - cellCoordinates[j].second)*(cellCoordinates[i].second - cellCoordinates[j].second) * dxC * dxC;
			rij = sqrt(rij);
			//if (abs(daYcorrArr[ida] - (y0 + icellCoordinates[jn].first * dyC)) < 0.5 * dyC && abs(daXcorrArr[ida] - (x0 + icellCoordinates[jn].second * dxC)) < 0.5 * dxC)
			if(rij < minDist)
			{
				minDist = rij;
				Hc_hgVector(ida) = jn;
				/*Hc_hgVector(im) = jn;    // *ns_statSize + stateIndex[io];  the index of the grid cell where there is observation
				im++;
				break;*/
			}
		}
	}

	return;
	//std::cout << std::endl<<" H Matrix: "<<std::endl << H_hMaxtirx << std::endl;
}
// read forcing data assimilation control file
void uebEnKFDA::readDaContr(const char* daconFile)   //, daControlV &inpDaControlV)
{
	std::ifstream pinFile(daconFile);
	char headerLine[256];
	pinFile.getline(headerLine, 256);   //skip header 
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d ", &es_enseSize);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d ", &es_pfSize_Fact);

	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d %d %d %d ", &ModelStartDateEns[0], &ModelStartDateEns[1], &ModelStartDateEns[2], &ModelStartDateEns[3]);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d %d %d %d ", &ModelEndDateEns[0], &ModelEndDateEns[1], &ModelEndDateEns[2], &ModelEndDateEns[3]);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d %d %d %d ", &forecastDateTimeEns[0], &forecastDateTimeEns[1], &forecastDateTimeEns[2], &forecastDateTimeEns[3]);

	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%s ", &xmrg1dEnsDir);
	//std::strcpy(xmrg1dEnsDir, "./outdaHF");
	char ueb_obsOutsideWS[256];
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%s ", &ueb_obsOutsideWS);
	if (strcmp(ueb_obsOutsideWS, "True") == 0 || strcmp(ueb_obsOutsideWS, "true") == 0)
		obsOutsideWS = true;
	else
		obsOutsideWS = false;
		
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &forcEnStdev);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &tempEnStdev);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &obsErrStdev);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &obsQErrStdev);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &daStatesStdev);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &daStatesStdev2);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d ", &qUpdateFreq);   //Q update frequency 
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &forcCorLength);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &dastateCorLength); 
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &tdecorrLength);      //temporal decorrelation length
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f %f ", &dyC, &dxC);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d ", &numObsPoints);   //10.6.17 for point site variables at obsr. (SNOTEL) stations 
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%d ", &daderivedType);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &dasnIndx);
	pinFile.getline(headerLine, 256, '\n');
	sscanf(headerLine, "%f ", &polyThreshold);
	/*std::cout << " " << es_enseSize << " " << es_pfSize_Fact << " " << forcEnStdev << " " << tempEnStdev << " " << obsErrStdev << " " << obsQErrStdev
	<< " " << daStatesStdev << " " << daStatesStdev2 << " " << forcCorLength << " " << dastateCorLength << " " << dyC << " " << " " << dxC
	<< " " << numObsPoints << " " << daderivedType << " " << dasnIndx << std::endl;*/
	for (int ic = 0; ic < 5; ic++) {
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%f ", &polyCoeff1[ic]);
		//daContArr[ida].polyCoeff[ic] = polyCoeff;
		//std::cout << polyCoeff1[ic] << "  ";
	}
	//std::cout << std::endl;
	for (int ic = 0; ic < 5; ic++) {
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%f ", &polyCoeff2[ic]);
		//daContArr[ida].polyCoeff[ic] = polyCoeff;
		//std::cout << polyCoeff2[ic] << "  ";
	}
	//std::cout << std::endl;	

	pinFile.close();
	return;
}
// read input text file and record datetime, and list of values; skip no data, get no data value from file
void uebEnKFDA::readTStextFile_multiVal(const char* inforcFile)       //, int &numdaPoints)  //, daControlV &inpDaControlV)
{
	FILE* inputFile = fopen(inforcFile, "r");
	if (!inputFile)
	{
		std::cout << "Error opening file: " << inforcFile << std::endl;
		return;
	}
	int nrecords = 0;
	int m_numObs;
	float noDataV = -9999;
	char commentLine[256];                    //string to read header line	
	fscanf(inputFile, "%f %d ", &noDataV, &m_numObs); //   get no data value, number of data cols
	if (m_numObs != numObsPoints)
	{
		std::cout << "Error the number of : " << m_numObs << " must be equal to number of obs points: " << numObsPoints << std::endl;
		std::getchar();
		return;
	}
	daYcorrArr.resize(m_numObs);
	daXcorrArr.resize(m_numObs);
	fgets(commentLine, 256, inputFile);   //skip remaining line
	for (int id = 0; id < m_numObs; id++)
		fscanf(inputFile, "%f %f ", &daYcorrArr[id], &daXcorrArr[id]); //   coordinates of data points
	fgets(commentLine, 256, inputFile);   //skip remaining contents of line

	fgets(commentLine, 256, inputFile);   //skip header line 
	int Year, Month, Day;
	double Hour, DTimeV;
	float Value;
	while (!feof(inputFile))
	{
		commentLine[0] = ' ';
		fgets(commentLine, 256, inputFile);
		if (commentLine[0] != ' ')
			++nrecords;
	}//while
	nRecs = nrecords;
	//strinpts = new inptimeseries[nrecords];                //assign memory to store data records
	daTcorrArr.resize(nrecords);
	daRegArray.resize(nrecords);
	for (int ir = 0; ir < nrecords; ir++)
		daRegArray[ir].resize(m_numObs);
	//
	rewind(inputFile);
	fgets(commentLine, 256, inputFile);   //no data value, number of data cols
	fgets(commentLine, 256, inputFile);   //coordinates of data points
	fgets(commentLine, 256, inputFile);   //skip header line 
	int inputRead = 0;
	for (int ir = 0; ir < nrecords; ir++)
	{
		fscanf(inputFile, "%d %d %d %lf %f ", &Year, &Month, &Day, &Hour, &Value);
		//std::cout << " hour " << Hour;
		if (fabs(Value - noDataV) > 0.1) {             //only copy data that is not no-data  ===>>> ***** 12.14.16; this needs revisioin 
			DTimeV = julian(Year, Month, Day, Hour);
			//std::cout << " hour julian " << std::setprecision(15)<< DTimeV;
			daTcorrArr[ir] = DTimeV;
			daRegArray[ir](0) = Value;
			for (int id = 1; id < m_numObs; id++)
				fscanf(inputFile, "%f ", &daRegArray[ir](id));
		}	 //
			 //fscanf(inputFile, " %*s\n");
		fgets(commentLine, 256, inputFile);   //skip remaining contents of line

		inputRead++;
	}

	std::cout << " Read file:  " << inforcFile << " number of lines: " << inputRead << std::endl;
	/*std::cout << " Read obs array at " << std::endl;
	for (int ir = 0; ir < nrecords; ir++)
	std::cout << " " << daTcorrArr[ir]<< " ";
	std::cout << std::endl;

	for (int ir = 0; ir < nrecords; ir++)
	{
	for (int id = 0; id < m_numObs; id++)
	std::cout << " " << daRegArray[ir][id] << " ";
	std::cout << std::endl;
	}
	std::getchar();*/

	fclose(inputFile);

	return;
}
// read input text file obs Q
void uebEnKFDA::readTStextFileQ(const char* inforcFile)       //, int &numdaPoints)  //, daControlV &inpDaControlV)
{
	FILE* inputFile = fopen(inforcFile, "r");
	if (!inputFile)
	{
		std::cout << "Error opening file: " << inforcFile << std::endl;
		return;
	}
	int nrecords = 0;
	float noDataV = -9999;
	char commentLine[256];                    //string to read header line	
	fscanf(inputFile, "%f ", &noDataV); //   get no data value, number of data cols
	fgets(commentLine, 256, inputFile);   //skip remaining line

	fgets(commentLine, 256, inputFile);   //skip header line 
	int Year, Month, Day;
	double Hour, DTimeV;
	float Value;
	while (!feof(inputFile))
	{
		commentLine[0] = ' ';
		fgets(commentLine, 256, inputFile);
		if (commentLine[0] != ' ')
			++nrecords;
	}//while
	 //strinpts = new inptimeseries[nrecords];                //assign memory to store data records
	daTcorrArrQ.resize(nrecords);
	daQstrArray.resize(nrecords);
	//
	rewind(inputFile);
	fgets(commentLine, 256, inputFile);   //no data value, number of data cols
	fgets(commentLine, 256, inputFile);   //skip header line 
	int inputRead = 0;
	for (int ir = 0; ir < nrecords; ir++)
	{
		fscanf(inputFile, "%d %d %d %lf %f \n", &Year, &Month, &Day, &Hour, &Value);
		//std::cout << " hour " << Hour;
		if (fabs(Value - noDataV) > 0.1) {             //only copy data that is not no-data  ===>>> ***** 12.14.16; this needs revisioin 
			DTimeV = julian(Year, Month, Day, Hour);
			//std::cout << " hour julian " << std::setprecision(15)<< DTimeV;
			daTcorrArrQ[ir] = DTimeV;
			daQstrArray(ir) = Value;
		}	 //
			 //fscanf(inputFile, " %*s\n");
			 //fgets(commentLine, 256, inputFile);   //skip remaining contents of line
		inputRead++;
	}
	std::cout << " Read file:  " << inforcFile << " number of lines: " << inputRead << std::endl;
	/*std::cout << " Read obs array at " << std::endl;
	for (int ir = 0; ir < nrecords; ir++)
	std::cout << " " << daTcorrArrQ[ir]<< " ";
	std::cout << std::endl;
	std::cout << " " << daQstrArray << " ";
	std::getchar();*/

	fclose(inputFile);

	return;
}
/*__host__ __device__*/
void uebEnKFDA::runEnKFopt(int ido, Eigen::VectorXf Z_obs, Eigen::RowVectorXf Xh_obsState, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> ensObservationErr, Eigen::RowVectorXf stateOutputArr, Eigen::RowVectorXf& stateOutputUpdate)         // bool NormalDist)   //float* &ensembleUpdateArr,
{
	float XhensMeanArr = Xh_obsState.mean();
	//ensemble anomaly
	Eigen::RowVectorXf Xh_ensAnomalyArr(es_enseSize);
	for (int ie = 0; ie < es_enseSize; ie++)
		Xh_ensAnomalyArr(ie) = Xh_obsState(ie) - XhensMeanArr;
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " obs grid: "<< ido <<" States in Obs space: " << std::endl;
		std::cout << Xh_obsState << " " << std::endl;
		std::cout << std::endl << " Ensemble mean in Obs space: " << std::endl;
		std::cout << XhensMeanArr << " " << std::endl;
		std::cout << std::endl << " Xh' Ensemble anomaly in Obs Space : " << std::endl;
		std::cout << Xh_ensAnomalyArr << " " << std::endl;
	}
	float Pzz_obsStateCov = Xh_ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pzz_obsStateCov /= (es_enseSize - 1);
	Pzz_obsStateCov += R_obsErrCov(ido, ido);		// .cast<float>();
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " R observation Error covariance matrix: " << std::endl;
		std::cout << R_obsErrCov(ido, ido) << " ";
	}
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pzz Innovation covariance matrix(Observation uncertainty) : " << std::endl;
		std::cout << Pzz_obsStateCov << std::endl;
	}
	//Pzz_i inverse of Observation uncertainty (innovation?) matrix
	float Pzz_i;
	if (Pzz_obsStateCov == 0) {
		std::cout << std::endl << " Warnning: div by zero! Press 'Enter' to continue with Kalman Gain = 0" << Pzz_obsStateCov << std::endl;
		std::ofstream debugOutputFile;
		debugOutputFile.open("debugOutput.txt", std::ios::app);
		debugOutputFile << " Warnning: div by zero !" << Pzz_obsStateCov << std::endl;
		debugOutputFile.close();
		std::getchar();
	}
	else
		Pzz_i = 1.0 / Pzz_obsStateCov;

	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pzzi Inverse of Pzz (Observation uncertainty / innovation covariance) matrix: " << std::endl;
		std::cout << Pzz_i << std::endl;
	}
	//
	float ensMeanArr = stateOutputArr.mean();
	//ensemble anomaly
	Eigen::RowVectorXf ensAnomalyArr(es_enseSize);
	for (int ie = 0; ie < es_enseSize; ie++)
		ensAnomalyArr(ie) = stateOutputArr(ie) - ensMeanArr;
	/*std::cout << std::endl << " States matrix: " << std::endl << stateOutputArr << " " << std::endl;
	std::cout << std::endl << " observed array: " << std::endl << Z_obs << " " << std::endl;
	*/
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " States matrix: " << std::endl;
		std::cout << stateOutputArr << " " << std::endl;
		std::cout << std::endl << " Ensemble mean: " << std::endl;
		std::cout << ensMeanArr << " " << std::endl;
		std::cout << std::endl << " Ensemble anomaly: " << std::endl;
		std::cout << ensAnomalyArr << " " << std::endl;
	}
	// Pxz state obs cross-covariance
	float Pxz_stateObsXCov = ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pxz_stateObsXCov /= (es_enseSize - 1);
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pxz state-obs cross-covariance matrix: " << std::endl;
		std::cout << Pxz_stateObsXCov << std::endl;
		std::cout << std::endl << " observed array: " << std::endl;
		std::cout << Z_obs << " " << std::endl;
		//
		std::cout << std::endl << "observation error: " << std::endl;
		std::cout << ensObservationErr << std::endl;
	}
	//residual y = z - HXb and R' = function of Yobs
	Eigen::RowVectorXf y_obsStateRresidual( es_enseSize);
	for (int ie = 0; ie < es_enseSize; ie++)
		y_obsStateRresidual(ie) = Z_obs(ido) - Xh_obsState(ie); //3.2.18 +ensObservationErr.cast<float>()(id, ie);  //added 10.13.17   
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " correction / residulal / Innovation (observed - model forecasted observation) y " << std::endl;
		std::cout << y_obsStateRresidual << std::endl;
	}
	//Kalman gain
	float K_kalGain;
	if (Pzz_obsStateCov == 0) {
		std::cout << std::endl << " Warnning: div by zero! Press 'enter' to continue with Kalman Gain = 0" << Pzz_obsStateCov << std::endl;
		K_kalGain = 0.0;                              //8.28.16 effect is unknown---so don't use update
		std::getchar();
		goto label1;
	}
	//K Kalman gain K = Pxz * Pzz^-1 = Pxz*Pb_R_1
	K_kalGain = Pxz_stateObsXCov * Pzz_i;
	//8.28.16 ----what if Kalman gain is = 0 for long time?
label1:
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " K Kalman gain: " << K_kalGain << std::endl;
	}
	//update state Xa = Xb + K*y, y = z-HXb
	Eigen::RowVectorXf X_state(es_enseSize);
	X_state = stateOutputArr + (K_kalGain * y_obsStateRresidual);
	//background Pxx = (1 / (es_enseSize - 1)) *ensAnomalyArr * ensAnomalyArr.transpose();
	float Pxx_stateCov = ensAnomalyArr * ensAnomalyArr.transpose();
	Pxx_stateCov /= (es_enseSize - 1);
	P_stateCov = Pxx_stateCov - K_kalGain * Pzz_obsStateCov * K_kalGain;
	//std::cout << std::endl << " Updated state matrix: " << std::endl << X_state << std::endl;
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Updated state matrix: " << std::endl;
		std::cout << X_state << std::endl;
		std::cout << std::endl << " Pxx Model background error covariance: " << std::endl;
		std::cout << Pxx_stateCov << std::endl;
		std::cout << std::endl << " Pxx_u Updated model error covariance: " << std::endl;
		std::cout << P_stateCov << std::endl;
	}
	//if SWE <= 0: 
	for (int ie = 0; ie < es_enseSize; ie++)
	{
		if (X_state(ie) < 0)
			X_state(ie) = 0.0;
		stateOutputUpdate(ie) = X_state(ie);
	}
	stateOutputUpdate(es_enseSize) = X_state.mean();

	if (uebDAsacrutpix7_debugout1 == 1)
	{
		std::cout << std::endl << " Updated state matrix mean: " << std::endl;
		std::cout << stateOutputUpdate << std::endl;
	}

	return;
}
//11.23.18 each obs points da without its observation (Leave one approach)
/*__host__ __device__*/
void uebEnKFDA::runEnKF_LeaveOneOut(int iog, Eigen::VectorXf Z_obs, std::vector<Eigen::RowVectorXf> Xh_obsState, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> ensObservationErr, Eigen::RowVectorXf stateOutputArr, Eigen::RowVectorXf& stateOutputUpdate)         // bool NormalDist)   //float* &ensembleUpdateArr,
{
	VectorXf XhensMeanArr(numObsPoints);
	for (int io = 0; io < numObsPoints; io++)
		XhensMeanArr(io) = Xh_obsState[io].mean();
	//ensemble anomaly
	Matrix<float, Dynamic, Dynamic, RowMajor> Xh_ensAnomalyArr(numObsPoints-1, es_enseSize);
	int ido = 0;
	for (int id = 0; id < numObsPoints; id++)
		if (id != iog)   //11.23.18 skip the grid cell iog
		{
			for (int ie = 0; ie < es_enseSize; ie++)
				Xh_ensAnomalyArr(ido, ie) = Xh_obsState[id](ie) - XhensMeanArr(id);
			ido++;
		}
	if (uebDAsacrutpix7_debugout2 == 1)	{
		std::cout << std::endl << " States in Obs space: " << std::endl;
		for (int id = 0; id < numObsPoints; id++)
			std::cout << Xh_obsState[id] << " " << std::endl;
		std::cout << std::endl << " Ensemble mean in Obs space: " << std::endl;
		std::cout << XhensMeanArr << " " << std::endl;
		std::cout << std::endl << " Xh' Ensemble anomaly in Obs Space : " << std::endl;
		std::cout << Xh_ensAnomalyArr << " " << std::endl;
	}
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>
		Pzz_obsStateCov = Xh_ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pzz_obsStateCov /= (es_enseSize - 1);
	ido = 0;
	int iro = 0;
	for (int id = 0; id < numObsPoints; id++)
	{
		iro = 0;
		if (id != iog)   //11.23.18 skip the grid cell iog
		{
			for (int ir = 0; ir < numObsPoints; ir++)
				if (ir != iog)   //11.23.18 skip the grid cell iog
				{
					Pzz_obsStateCov(ido, iro) += R_obsErrCov(id, ir);		// .cast<float>();
					iro++;
				}
			ido++;
		}
	}
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " R observation Error covariance matrix: " << std::endl;
		std::cout << R_obsErrCov << " ";
	}
	if (uebDAsacrutpix7_debugout2 == 1)	{
		std::cout << std::endl << " Pzz Innovation covariance matrix(Observation uncertainty) : " << std::endl;
		std::cout << Pzz_obsStateCov << std::endl;
	}
	//Pzz_i inverse of Observation uncertainty (innovation?) matrix
	Matrix<float, Dynamic, Dynamic, RowMajor> Pzz_i(numObsPoints-1, numObsPoints-1);
	float Pzz_det = Pzz_obsStateCov.determinant();
	if (Pzz_det == 0) {
		std::cout << std::endl << " Warnning: zero determinant of matrix! Press 'Enter' to continue with Kalman Gain = 0" << Pzz_det << std::endl;
		std::ofstream debugOutputFile;
		debugOutputFile.open("debugOutput.txt", std::ios::app);
		debugOutputFile << " Warnning: zero determinant of matrix!" << Pzz_det << std::endl;
		debugOutputFile.close();
		std::getchar();
	}
	else
		Pzz_i = Pzz_obsStateCov.inverse();

	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pzzi Inverse of Pzz (Observation uncertainty / innovation covariance) matrix: " << std::endl;
		std::cout << Pzz_i << std::endl;
	}
	//
	float ensMeanArr = stateOutputArr.mean();
	//ensemble anomaly
	Eigen::RowVectorXf ensAnomalyArr(es_enseSize);
	for (int ie = 0; ie < es_enseSize; ie++)
		ensAnomalyArr(ie) = stateOutputArr(ie) - ensMeanArr;
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " States matrix: " << std::endl;
		std::cout << stateOutputArr << " " << std::endl;
		std::cout << std::endl << " Ensemble mean: " << std::endl;
		std::cout << ensMeanArr << " " << std::endl;
		std::cout << std::endl << " Ensemble anomaly: " << std::endl;
		std::cout << ensAnomalyArr << " " << std::endl;
	}
	// Pxz state obs cross-covariance
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>
		Pxz_stateObsXCov = ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pxz_stateObsXCov /= (es_enseSize - 1);
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pxz state-obs cross-covariance matrix: " << std::endl;
		std::cout << Pxz_stateObsXCov << std::endl;
		std::cout << std::endl << " observed array: " << std::endl;
		std::cout << Z_obs << " " << std::endl;
		//
		std::cout << std::endl << "observation error: " << std::endl;
		std::cout << ensObservationErr << std::endl;
	}
	//residual y = z - HXb and R' = function of Yobs
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> y_obsStateRresidual(numObsPoints-1, es_enseSize);
	ido = 0;
	for (int id = 0; id < numObsPoints; id++)
		if (id != iog)   //11.23.18 skip the grid cell iog
		{
			for (int ie = 0; ie < es_enseSize; ie++)
				y_obsStateRresidual(ido, ie) = Z_obs(id) - Xh_obsState[id](ie) + ensObservationErr(id, ie);  //added 10.13.17   
			ido++;
		}
	if (uebDAsacrutpix7_debugout2 == 1)	{
		std::cout << std::endl << " correction / residulal / Innovation (observed - model forecasted observation) y " << std::endl;
		std::cout << y_obsStateRresidual << std::endl;
	}
	//Kalman gain
	Eigen::RowVectorXf K_kalGain(numObsPoints-1);
	if (Pzz_det == 0) {
		std::cout << std::endl << " Warnning: zero determinant of matrix! Press 'enter' to continue with Kalman Gain = 0" << Pzz_det << std::endl;
		K_kalGain.setZero();                                 //8.28.16 effect is unknown---so don't use update
		std::getchar();
		goto label1;
	}
	//K Kalman gain K = Pxz * Pzz^-1 = Pxz*Pb_R_1
	K_kalGain = Pxz_stateObsXCov * Pzz_i;
	//8.28.16 ----what if Kalman gain is = 0 for long time?
label1:
	if (uebDAsacrutpix7_debugout2 == 1)	{
		std::cout << std::endl << " K Kalman gain: " << std::endl;
		std::cout << K_kalGain << std::endl;
	}
	//update state Xa = Xb + K*y, y = z-HXb
	Eigen::RowVectorXf X_state(es_enseSize);
	X_state = stateOutputArr + (K_kalGain * y_obsStateRresidual);
	//background Pxx = (1 / (es_enseSize - 1)) *ensAnomalyArr * ensAnomalyArr.transpose();
	float Pxx_stateCov = ensAnomalyArr * ensAnomalyArr.transpose();
	Pxx_stateCov /= (es_enseSize - 1);
	P_stateCov = Pxx_stateCov - K_kalGain * Pzz_obsStateCov * K_kalGain.transpose();
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Updated state matrix: " << std::endl;
		std::cout << X_state << std::endl;
		std::cout << std::endl << " Pxx Model background error covariance: " << std::endl;
		std::cout << Pxx_stateCov << std::endl;
		std::cout << std::endl << " Pxx_u Updated model error covariance: " << std::endl;
		std::cout << P_stateCov << std::endl;
	}
	//if SWE <= 0: 
	for (int ie = 0; ie < es_enseSize; ie++)
	{
		if (X_state(ie) < 0)
			X_state(ie) = 0.0;
		stateOutputUpdate(ie) = X_state(ie);
	}
	stateOutputUpdate(es_enseSize) = X_state.mean();

	if (uebDAsacrutpix7_debugout1 == 1)
	{
		std::cout << std::endl << " Updated state matrix mean: " << std::endl;
		std::cout << stateOutputUpdate << std::endl;
	}

	return;
}
/*__host__ __device__*/
void uebEnKFDA::runEnKF(Eigen::VectorXf Z_obs, std::vector<Eigen::RowVectorXf> Xh_obsState, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> ensObservationErr, Eigen::RowVectorXf stateOutputArr, Eigen::RowVectorXf& stateOutputUpdate)         // bool NormalDist)   //float* &ensembleUpdateArr,
{
	VectorXf XhensMeanArr(numObsPoints);
	for (int io = 0; io < numObsPoints; io++)
		XhensMeanArr(io) = Xh_obsState[io].mean();
	//ensemble anomaly
	Matrix<float, Dynamic, Dynamic, RowMajor> Xh_ensAnomalyArr(numObsPoints, es_enseSize);
	for (int id = 0; id < numObsPoints; id++)
		for (int ie = 0; ie < es_enseSize; ie++)
			Xh_ensAnomalyArr(id, ie) = Xh_obsState[id](ie) - XhensMeanArr(id);
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " States in Obs space: " << std::endl;
		for (int id = 0; id < numObsPoints; id++)
			std::cout << Xh_obsState[id] << " " << std::endl;
		std::cout << std::endl << " Ensemble mean in Obs space: " << std::endl;
		std::cout << XhensMeanArr << " " << std::endl;
		std::cout << std::endl << " Xh' Ensemble anomaly in Obs Space : " << std::endl;
		std::cout << Xh_ensAnomalyArr << " " << std::endl;
	}
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>
		Pzz_obsStateCov = Xh_ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pzz_obsStateCov /= (es_enseSize - 1);
	Pzz_obsStateCov += R_obsErrCov;		// .cast<float>();
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " R observation Error covariance matrix: " << std::endl;
		std::cout << R_obsErrCov << " ";
	}
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pzz Innovation covariance matrix(Observation uncertainty) : " << std::endl;
		std::cout << Pzz_obsStateCov << std::endl;
	}
	//Pzz_i inverse of Observation uncertainty (innovation?) matrix
	Matrix<float, Dynamic, Dynamic, RowMajor> Pzz_i(numObsPoints, numObsPoints);
	float Pzz_det = Pzz_obsStateCov.determinant();
	if (Pzz_det == 0) {
		std::cout << std::endl << " Warnning: zero determinant of matrix! Press 'Enter' to continue with Kalman Gain = 0" << Pzz_det << std::endl;
		std::ofstream debugOutputFile;
		debugOutputFile.open("debugOutput.txt", std::ios::app);
		debugOutputFile << " Warnning: zero determinant of matrix!" << Pzz_det << std::endl;
		debugOutputFile.close();
		std::getchar();
	}
	else
		Pzz_i = Pzz_obsStateCov.inverse();

	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pzzi Inverse of Pzz (Observation uncertainty / innovation covariance) matrix: " << std::endl;
		std::cout << Pzz_i << std::endl;
	}
	//
	float ensMeanArr = stateOutputArr.mean();
	//ensemble anomaly
	Eigen::RowVectorXf ensAnomalyArr(es_enseSize);
	for (int ie = 0; ie < es_enseSize; ie++)
		ensAnomalyArr(ie) = stateOutputArr(ie) - ensMeanArr;
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " States matrix: " << std::endl;
		std::cout << stateOutputArr << " " << std::endl;
		std::cout << std::endl << " Ensemble mean: " << std::endl;
		std::cout << ensMeanArr << " " << std::endl;
		std::cout << std::endl << " Ensemble anomaly: " << std::endl;
		std::cout << ensAnomalyArr << " " << std::endl;
	}
	// Pxz state obs cross-covariance
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>
		Pxz_stateObsXCov = ensAnomalyArr * Xh_ensAnomalyArr.transpose();
	Pxz_stateObsXCov /= (es_enseSize - 1);
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Pxz state-obs cross-covariance matrix: " << std::endl;
		std::cout << Pxz_stateObsXCov << std::endl;
		std::cout << std::endl << " observed array: " << std::endl;
		std::cout << Z_obs << " " << std::endl;
		//
		std::cout << std::endl << "observation error: " << std::endl;
		std::cout << ensObservationErr << std::endl;
	}
	//residual y = z - HXb and R' = function of Yobs
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> y_obsStateRresidual(numObsPoints, es_enseSize);
	for (int id = 0; id < numObsPoints; id++)
		for (int ie = 0; ie < es_enseSize; ie++)
			y_obsStateRresidual(id, ie) = Z_obs(id) - Xh_obsState[id](ie) + ensObservationErr(id, ie);  //added 10.....18 
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " correction / residulal / Innovation (observed - model forecasted observation) y " << std::endl;
		std::cout << y_obsStateRresidual << std::endl;
	}
	//Kalman gain
	Eigen::RowVectorXf K_kalGain(numObsPoints);
	if (Pzz_det == 0) {
		std::cout << std::endl << " Warnning: zero determinant of matrix! Press 'enter' to continue with Kalman Gain = 0" << Pzz_det << std::endl;
		K_kalGain.setZero();                                 //8.28.16 effect is unknown---so don't use update
		std::getchar();
		goto label1;
	}
	//K Kalman gain K = Pxz * Pzz^-1 = Pxz*Pb_R_1
	K_kalGain = Pxz_stateObsXCov * Pzz_i;
	//8.28.16 ----what if Kalman gain is = 0 for long time?
label1:
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " K Kalman gain: " << std::endl;
		std::cout << K_kalGain << std::endl;
	}
	//update state Xa = Xb + K*y, y = z-HXb
	Eigen::RowVectorXf X_state(es_enseSize);
	X_state = stateOutputArr + (K_kalGain * y_obsStateRresidual);
	//background Pxx = (1 / (es_enseSize - 1)) *ensAnomalyArr * ensAnomalyArr.transpose();
	float Pxx_stateCov = ensAnomalyArr * ensAnomalyArr.transpose();
	Pxx_stateCov /= (es_enseSize - 1);
	P_stateCov = Pxx_stateCov - K_kalGain * Pzz_obsStateCov * K_kalGain.transpose();
	if (uebDAsacrutpix7_debugout2 == 1)
	{
		std::cout << std::endl << " Updated state matrix: " << std::endl;
		std::cout << X_state << std::endl;
		std::cout << std::endl << " Pxx Model background error covariance: " << std::endl;
		std::cout << Pxx_stateCov << std::endl;
		std::cout << std::endl << " Pxx_u Updated model error covariance: " << std::endl;
		std::cout << P_stateCov << std::endl;
	}
	//if SWE <= 0: 
	for (int ie = 0; ie < es_enseSize; ie++)
	{
		if (X_state(ie) < 0)
			X_state(ie) = 0.0;
		stateOutputUpdate(ie) = X_state(ie);
	}
	stateOutputUpdate(es_enseSize) = X_state.mean();

	if (uebDAsacrutpix7_debugout1 == 1)
	{
		std::cout << std::endl << " Updated state matrix mean: " << std::endl;
		std::cout << stateOutputUpdate << std::endl;
	}

	return;
}
void uebEnKFDA::updateDaArr(int& startIndexDA)     //float** RegArray, double* tvarArr, daControlV inpDaControlV)          //, int nDataPoints) 
{
	//need to call each variable array as each array has to be copied separately to device array in cuda	
	int io = 0;
	float xcorVal, ycorVal;
	float expPower;
	for (int iot = 0; iot < numObsPoints; iot++)
	{
		xcorVal = daRegArray[startIndexDA][iot];
		if (daderivedType == 1)
		{

			if (xcorVal < polyThreshold)
			{
				ycorVal = polyCoeff1[0] * powf(xcorVal, 4) + polyCoeff1[1] * powf(xcorVal, 3) +
					polyCoeff1[2] * powf(xcorVal, 2) + polyCoeff1[3] * powf(xcorVal, 1) +
					polyCoeff1[4] * powf(xcorVal, 0);
			}
			else
			{
				ycorVal = polyCoeff2[0] * powf(xcorVal, 4) + polyCoeff2[1] * powf(xcorVal, 3) +
					polyCoeff2[2] * powf(xcorVal, 2) + polyCoeff2[3] * powf(xcorVal, 1) +
					polyCoeff2[4] * powf(xcorVal, 0);
			}
			Z_obs[iot] = dasnIndx *  ycorVal;
			//std::cout << "obs at SNOTEL: " << Z_obs[iot] << std::endl;
		}
		else
			Z_obs[iot] = xcorVal;
	}
	daTime = daTcorrArr[startIndexDA];
	startIndexDA = startIndexDA + 1;
	//std::cout << " updated DA arr at : " << startIndexDA<<" time "<< daTime << std::endl;
	//std::cout << "obs at SNOTEL: " << std::endl;
	//std::cout << Z_obs << std::endl;

	return;
}
//12.16.17 for da out.........//output control file: details of netcdf outputs 
void uebEnKFDA::readDAOutputControl(const char* outputconFile)      //, std::vector<ncdaOutput> &daOut, std::vector<ncdaOutput> &ensOut)
{
	std::ifstream poutFile(outputconFile);
	char headerLine[256];
	int nout;
	//istringstream valueLine;	
	poutFile.getline(headerLine, 256);   //skip header	
										 //data assimilation output
	poutFile.getline(headerLine, 256);
	sscanf(headerLine, "%d ", &nout);
	//ndaout = nout;
	//daOut = new ncOutput[nout];
	daOutArr.resize(nout);
	for (int i = 0; i < nout; i++)
	{
		poutFile.getline(headerLine, 256);
		sscanf(headerLine, "%s %s %s", &daOutArr[i].symbol, &daOutArr[i].outfName, &daOutArr[i].units);
		for (int sindx = 0; sindx < 9; sindx++)     //TODO: try to do without looping?  
		{
			if (strcmp(daOutArr[i].symbol, uebDAsacrutpix7outStates[sindx]) == 0)
			{
				daOutArr[i].outvarIndex = sindx;
				break;
			}
		}
	}
	//ensemble output without data assimilastion
	poutFile.getline(headerLine, 256);
	sscanf(headerLine, "%d ", &nout);
	//nensout = nout;
	//ensOut = new ncOutput[nout];
	daEnsArr.resize(nout);
	for (int i = 0; i < nout; i++)
	{
		poutFile.getline(headerLine, 256);
		sscanf(headerLine, "%s %s %s", &daEnsArr[i].symbol, &daEnsArr[i].outfName, &daEnsArr[i].units);
		for (int sindx = 0; sindx < 9; sindx++)     //TODO: try to do without looping?  
		{
			if (strcmp(daEnsArr[i].symbol, uebDAsacrutpix7outStates[sindx]) == 0)
			{
				daEnsArr[i].outvarIndex = sindx;
				break;
			}
		}
	}
	//
	poutFile.close();
	//paoutFile.close();
	return;
}
//***************************** JULIAN () ****************************
//             To convert the real date to julian date
// YJS The Julian are change to a new version to take the Leap Yean into consideration
//    in the old version, there are 365 days each year.
//     FUNCTION JULIAN(MONTH,DAY)
/*__host__ __device__*/
int uebEnKFDA::julian(int yy, int mm, int dd)
{
	int julian;
	int mmstrt[12] = { 0,31,59,90,120,151,181,212,243,273,304,334 };
	int jday = mmstrt[mm - 1] + dd;
	int ileap = yy - ((int)(yy / 4)) * 4;
	if ((ileap == 0) && (mm >= 3))
		jday = jday + 1;
	julian = jday;
	return julian;
}
//The following were copied from functions.f90 # 6.8.13
//THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT time VALUE
//CAN BE IN ANY UT-LIKE time SCALE (UTC, UT1, TT, ETC.) - OUTPUT. //JULIAN DATE WILL HAVE SAME BASIS.  
//ALGORITHM BY FLIEGEL AND //VAN FLANDERN. //SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
//I = YEAR (IN) //M = MONTH NUMBER (IN) //K = DAY OF MONTH (IN) //H = UT HOURS (IN) //TJD = JULIAN DATE (OUT)
/*__host__ __device__*/
double uebEnKFDA::julian(int I, int M, int K, double H)
{
	double TJD, JD;
	//JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
	JD = K - 32075 + 1461 * (I + 4800 + (M - 14) / 12) / 4 + 367 * (M - 2 - (M - 14) / 12 * 12) / 12 - 3 * ((I + 4900 + (M - 14) / 12) / 100) / 4;
	TJD = JD - 0.5 + H / 24.0;
	//##%^_TBC 6.8.13 //powf(10,0) in place of D0
	return TJD;
}

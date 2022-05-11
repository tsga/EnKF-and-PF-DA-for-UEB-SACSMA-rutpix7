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

#include "uebDAsrPFuebpgdecls.h"
 
uebCellDA::uebCellDA(float paramArr[32], int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT,
	double  UTCoffset, int inpDailyorSubd, int oStride)
{
	float siteVarInitcondefaults[32] = { 0.0, 0.0, 0.0, 0.0, 1.0, 100000.0, 0.1, 0.0, 0.0,
		0.0, 6.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.98, 5.712903,
		4.350000, 6.890322, 8.660001, 8.938710, 10.010000, 9.541936,
		9.038710, 7.160001, 8.106450, 5.923332, 5.058064, -9999.0, 111.00 };
	//float paramArr[32];
	setConstantValues();
	//set parameters	
	//readParams(inpFile, paramArr);
	/*std::cout<<"param read..\n ");
	for(int i=0;i<npar;i++)
	std::cout<<"%f ",parvalArray[i]);	*/
	setParams(paramArr);	
	//model run settings
	setModelRun_Settings(startDate, endDate, startHour, endHour, modeldT, UTCoffset, inpDailyorSubd, oStride);
	//site vars and intitial conditions 
	setSiteVars_and_Initconds(siteVarInitcondefaults);
        //std::cout<<"UEBCell initialized"<<std::endl;	
}

uebCellDA::uebCellDA()
{
	float paramDefaults[32] = {0, 0, 3, -1, 0.98, 2.09, 2, 0.01, 337, 1700, 0.05, 20, 0.1, 0.85,
		0.65, 0.278, 1.11, 0.0654, 1, 0.001, 0.98, 0.5, 0, 0.5, 0.004626286,
		0.25, 0.5, 0.857143, 0.16, 0.5, 0.8, 2.4 };
	float siteVarInitcondefaults[32] = { 0.0, 0.0, 0.0, 0.0, 1.0, 100000.0, 0.1, 0.0, 0.0,
		0.0, 6.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.98, 5.712903,
		4.350000, 6.890322, 8.660001, 8.938710, 10.010000, 9.541936,
		9.038710, 7.160001, 8.106450, 5.923332, 5.058064, -9999.0, 111.00 };
	setConstantValues();
	setParams(paramDefaults);	
	//model run settings
	//Default run settings 
	int startDate[3] = { 2009, 10, 1 }, endDate[3] = { 2010, 6, 1 };
	double startHour = 0.0, endHour = 0.0, modeldT = 1.0, UTCoffset = -7;
	int inpDailyorSubd = 0, outStride = 4;     //0 subdaily input
	setModelRun_Settings(startDate, endDate, startHour, endHour, modeldT, UTCoffset, inpDailyorSubd, outStride);
	//site vars and intitial conditions 
	setSiteVars_and_Initconds(siteVarInitcondefaults);
	//std::cout<<"UEBCell initialized"<<std::endl;	
}
/*
uebCellDA::uebCellDA(uebCellDA& uCell0)
{
	setConstantValues();
	setParams(uCell0.paramValues);
	//site vars and intitial conditions 
	setSiteVars_and_Initconds(uCell0.statesiteValues);
	//model run settings	
	setModelRun_Settings(uCell0.modelStartDate, uCell0.modelEndDate, uCell0.modelStartHour, uCell0.modelEndHour, uCell0.modelDT, uCell0.UTCOffset, uCell0.inpDailyorSubdaily, uCell0.outtStride);
	nstepinaDay = uCell0.nstepinaDay;
	if (uCell0.tsprevday)
	{
		tsprevday = new float[nstepinaDay];
		//  Initialize Tsbackup and TaveBackup
		for (int i = 0; i < nstepinaDay; i++)
		{
			tsprevday[i] = uCell0.tsprevday[i];
		}
	}
	else tsprevday = NULL;
	if (uCell0.taveprevday)
	{
		taveprevday = new float[nstepinaDay];
		for (int i = 0; i< nstepinaDay; i++)
		{
			taveprevday[i] = uCell0.taveprevday[i];
		}
	}
	else taveprevday = NULL;

}

uebCellDA& uebCellDA::operator= (uebCellDA& uCell0)
{
	if (this != &uCell0)
	{   
		setConstantValues();
		setParams(uCell0.paramValues);
		//site vars and intitial conditions 
		setSiteVars_and_Initconds(uCell0.statesiteValues);
		//model run settings	
		setModelRun_Settings(uCell0.modelStartDate, uCell0.modelEndDate, uCell0.modelStartHour, uCell0.modelEndHour, uCell0.modelDT, uCell0.UTCOffset, uCell0.inpDailyorSubdaily, uCell0.outtStride);
		nstepinaDay = uCell0.nstepinaDay;
		delete[] tsprevday;
		delete[] taveprevday;
		if (uCell0.tsprevday)
		{
			tsprevday = new float[nstepinaDay];
			//  Initialize Tsbackup and TaveBackup
			for (int i = 0; i < nstepinaDay; i++)
			{
				tsprevday[i] = uCell0.tsprevday[i];
			}
		}
		else tsprevday = NULL;
		if (uCell0.taveprevday)
		{
			taveprevday = new float[nstepinaDay];
			for (int i = 0; i< nstepinaDay; i++)
			{
				taveprevday[i] = uCell0.taveprevday[i];
			}
		}
		else taveprevday = NULL;
	}
	
	return *this;
}*/

uebCellDA::~uebCellDA()
{
	/*delete []tsprevday;
	delete []taveprevday;
	tsprevday = NULL;
	taveprevday = NULL;*/
	/*for(int i=0; i<71;i++)
		delete[] OutVarValues[i];
	delete []OutVarValues;*/
}

__host__ __device__ 
void uebCellDA::setConstantValues()
{
	//defalut is accumulation zone is false
	accumulationZone = false;
	//initialize 
	T_0 = 0.0;				// Temperature of freezing (0 C)
	T_k = 273.15;			// Temperature to convert C to K (273.15)
	SB_c = 2.041334e-7;      // Stefan boltzman constant (2.041334e-7 KJ/m^2-hr-K^4) #corrected 12.23.14
	H_f = 333.5;			// Heat of fusion (333.5 KJ= kg)
	Hne_u = 2834.0;		// Heat of Vaporization (Ice to Vapor, 2834 KJ= kg)
	C_w = 4.18;			// Water Heat Capacity (4.18 KJ/ kg-C)
	C_s = 2.09;			// Ice heat capacity (2.09 KJ/ kg-C)
	C_p = 1.005;			// Air Heat Capacity (1.005 KJ= kg= K)
	Ra_g = 287.0;			// Ideal Gas constant for dry air (287 J= kg= K)
	K_vc = 0.4;			// Von Karmans constant (0.4)
	Hs_f = 3600.0;		// Factor to convert = s into = hr (3600)
	Rho_i = 917.0;			// Density of Ice (917 kg= m^3)
	Rho_w = 1000.0;		// Density of Water (1000 kg= m^3)
	Gra_v = 9.81;			// Gravitational acceleration (9.81 m= s^2)
	W1da_y = 0.261799;		// Daily frequency (2pi= 24 hr 0.261799 radians= hr) 
	Io = 4914.0;            //  Solar constant  Kj/m^2/hr
	//pi copied from snowdxv.f90
	P_i = 3.141592653589793238462643383279502884197169399375105820974944592308;		// Pi
	//data for pred-corr
	wtol = 0.025;
	utol = 2000.0;
	//from TURBFLUX()
	tol = 0.001;
	nitermax = 20;
	ncitermax = 21;
	// flag to write warnings,...etc
	snowdgtvariteflag = 0;
	snowdgtvariteflag2 = 0;  // 0;
	snowdgtvariteflag3 = 0;
	snowdgt_outflag = 0;
	radwarnflag = 0;

	//added 9.16.13
	iflag[4] = 4;
	
	//inpDailyorSubdaily = 0;       // 0: values given at each (sub-daily time steps); 1: daily values
	uebCellDAX = 0;
	uebCellDAY = 0;
	Tsk_save = 273.16, Tssk_old = 273.16, Tsavek_old = 273.16, Tsavek_ave = 273.16, Tssk_ave = 273.16/*added 6.7.13*/;

	//## these were copied from snowdgtv, not clear where they are being used 
	/*fStab = -9999;
	Tref = -9999;
	iTsMethod = 4;
	//#_This is not clear 8.28.13
	windfl = 0;	*/
	return;
}

// functions to read params
void uebCellDA::readParams(const char* inpFile, float Params[32])
{
	std::ifstream pinFile(inpFile);
	char headerLine[256];	
	pinFile.getline(headerLine, 256, '\n'); //skip header
	for (int i = 0; i < 32; i++)
	{
		pinFile.getline(headerLine, 256, '\n');
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%f ", &Params[i]);
	}
	pinFile.close();
	return;
}

__host__ __device__ 
void uebCellDA::setParams(float Params[32])
{
	//copy params class variables
	for (int i = 0; i < 32; i++)
		paramValues[i] = Params[i];
//5.2.15 from snowdv
	//  Mapping from parameters read to UEB internal interpretation which follows UEBVeg scheme
	irad = (int)paramValues[0];
	ireadalb = (int)paramValues[1];
	for (int i = 0; i<11; i++)
		Param[i] = paramValues[i + 2];
	for (int i = 12; i<18; i++)
		Param[i] = paramValues[i + 1];
	Param[18] = -9999;
	Param[19] = -9999;
	Param[20] = paramValues[19];
	for (int i = 22; i<32; i++)
		Param[i] = paramValues[i - 2];
	bca = paramValues[30];
	bcc = paramValues[31];
	return;
}
//copy site variables and state intitial conditions at a grid (ueb cell)
__host__ __device__ 
void uebCellDA::setSiteVars_and_Initconds(float SiteVars[32])
{
	//
	modisAlbedoFact = 0.0;   // albedo multiplier for modis

	for (int i = 0; i < 32; i++)
		statesiteValues[i] = SiteVars[i];
	//copy initial conditions
//==============================================changes for new conf 5.1.15
	//from snowdv
	//copied from paramsiteinitial
	for (int i = 0; i < 4; i++)
		statev[i] = statesiteValues[i];
	sitev[0] = statesiteValues[4];
	sitev[1] = statesiteValues[5];
	for (int i = 3; i < 9; i++)
		sitev[i] = statesiteValues[i + 3];
	slope = statesiteValues[12];
	azi = statesiteValues[13];
	lat = statesiteValues[14];

	Param[11] = statesiteValues[15];
	//subalb=statesiteValues[15]
	sitev[9] = statesiteValues[16];
	subtype = (int)statesiteValues[16];
	Param[21] = statesiteValues[17];
	//gsurf = statesiteValues[17]
	for (int i = 0; i < 12; i++)
		dtbar[i] = statesiteValues[i + 18];
	ts_last = statesiteValues[30];
	lon = statesiteValues[31];

	if (subtype == 0 || subtype == 3)
		WGT = 0.0;
	else
		WGT = 1.0;

	if (subtype != 3)  //  Only do this work for non accumulation cells where model is run
	{
		//  Initialize Tsbackup and TaveBackup
		for (int i = 0; i < nstepinaDay; i++)
		{
			tsprevday[i] = -9999.0;
			taveprevday[i] = -9999.0;
		}
		// Take surface temperature as 0 where it is unknown the previous time step
		// This is for first day of the model to get the force restore going
		//#$#$#$#$#_is this all the time steps or the last time?
		if (ts_last <= -9999)
			//for(int i =0;i< nstepinaDay;i++)
			tsprevday[nstepinaDay - 1] = 0;
		else
			//for(int i =0;i< nstepinaDay;i++)
			tsprevday[nstepinaDay - 1] = ts_last;
		// compute Ave.Temp for previous day 
		Us = statev[0];                                  // Ub in UEB
		Ws = statev[1];                                  // W in UEB
		Wc = statev[3];                                  // Canopy SWE
		Apr = sitev[1];                                   // Atm. Pressure  [PR in UEB]
		cg = Param[3];                                   // Ground heat capacity [nominally 2.09 KJ/kg/C]
		rhog = Param[7];                                   // Soil Density [nominally 1700 kg/m^3]
		de = Param[10];                                  // Thermally active depth of soil (0.1 m)

		//this are for coudiness computation
		//6.10.13
		as = Param[27];
		bs = Param[28];

		tave = TAVG(Us, Ws + WGT, Rho_w, C_s, T_0, rhog, de, cg, H_f);
		//for(int i =0;i< nstepinaDay;i++)
		taveprevday[nstepinaDay - 1] = tave;
		//  initialize variables for mass balance
		Ws1 = statev[1];
		Wc1 = statev[3];
		cumP = 0.0;
		cumEs = 0.0;
		cumEc = 0.0;
		cumMr = 0.0;
		cumGm = 0.0;
		cumEg = 0.0;
		dStorage = 0.0;
		errMB = 0.0;
	} //  end the skip block done only for accumulation cells	
	Tmin = 0.0;
	Tmax = 0.0;

	return;
}
__host__ __device__ 
void uebCellDA::setModelRun_Settings(int startDate[3], int endDate[3], double  startHour, double  endHour, double  modeldT, double  UTCoffset, int inpDailyorSubd, int oStride)
{
	//for EnKF
	startIndexDA = 0;
	startIndexDAQ = 0;
	forcEnStdev = 0.1;    //forcing ensemble standard deviation except temperature
	tempEnStdev = 1.0;
	tdecorLength = 24.0;
	corrFact1 = 0.5;	// 1.0 - modelDT / tdecorLength;
	corrFact2 = 0.5;	// sqrtf(1.0 - corrFact1 * corrFact1);    //1.0 - corrFact1;  // 

	daAssimlate = true;
	updateDaArray = true;
	updateDaArrayQ = true;

	for (int i = 0; i<3; i++)
	{
		modelStartDate[i] = startDate[i];
		modelEndDate[i] = endDate[i];
	}
	modelStartHour = startHour;
	modelEndHour = endHour;
	modelDT = modeldT;	
	UTCOffset = UTCoffset;

	//inpDailyorSubdaily = inpDailyorSubd;
	
//5.2.15 from snowdv
	//  FIXME: what if the result is fractional
	//  time steps must divide exactly in to a day because we use logic that requires the values from the same time
	//  step on the previous day.  Consider in future making the specification of time step as number of time
	//  steps in a day, not modeldt to ensure this modeldt is recalculated based on the int timesteps in a day
	//  assumption: number of model timesteps in a day must be an int                    
	nstepinaDay = (int)(24.0 / modelDT + 0.5);  // closest rounding
	modelDT = 24.0 / nstepinaDay;
	//tsprevday = new float[nstepinaDay];
	//taveprevday = new float[nstepinaDay]; 

	// Variables to keep track of which time step we are in and which netcdf output file we are in

	istep = 0;  // time step initiated as 0
	// map on to old UEB names
	Year = modelStartDate[0];     //7.16.16
	Month = modelStartDate[1];
	Day = modelStartDate[2];
	sHour = modelStartHour;
	currentModelDateTime = julian(Year, Month, Day, sHour);
	//double dlastD = Day, dlastH = dHour;
	modelSpan = julian(modelEndDate[0], modelEndDate[1], modelEndDate[2], modelEndHour) - julian(modelStartDate[0], modelStartDate[1], modelStartDate[2], modelStartHour); //no of days in model span
	//model time steps																																									   //model time steps
	numTotalTimeSteps = (int)ceil(modelSpan*(24 / modelDT));

	// for block-time simulation (for GPU)
	numSimTimeSteps = 24;
	/*if (inpDailyorSubdaily == 0)
		numSimTimeSteps = 24;                         //number of time steps in one simulation batch, i.e. inputs are read for this number of time steps to reduce repeated disk access
	else
		numSimTimeSteps = 24 * nstepinaDay;*/

	for (int i = 0; i < 13; i++)
	{
		startIndex[i] = 0;
		ncReadStart[i] = 0;
	}
	//tEnd = 0;
	outtStride = oStride;
	timeSeriesIndex = 0;                  //this changes to 1 when a forcing that is applicable for the whole model is read - --- forcing time series from text file need to be read only once
	//allocate memory for output array	
	//OutVarValues = new float *[71];
	//for (int i = 0; i < 71; i++)
	//OutVarValues = new float[71*numTotalTimeSteps];         // *outtStride];
	//std::cout << "number of t " << numTotalTimeSteps << std::endl;
	/*tsprevday.clear();
	tsprevday.resize(nstepinaDay);
	taveprevday.clear();
	taveprevday.resize(nstepinaDay);
	for (int i = 0; i < nstepinaDay; i++)
	{
		tsprevday[i] = -9999.0;		
		taveprevday[i] = -9999.0;
	}	*/
// 5.2.15 from snowdv
	// calculating model end date-time in julian date
	dHour = modelEndHour;
	EJD = julian(modelEndDate[0], modelEndDate[1], modelEndDate[2], dHour);

	return;
}
//7.22.16
__host__ __device__
void uebCellDA::setInitialEnsembleStates(int nEns)   //, const char* forcName)
{
	corrFact1 = fabs(1.0 - modelDT / tdecorLength);   //avoid very small -ve 
	corrFact2 = sqrtf(1.0 - corrFact1 * corrFact1);    //1.0 - corrFact1;  // 

	if (nEns > 59)
	{
		std::cout << " Error! Ensemble size bigger than 60 is not supported." << std::endl;
		std::getchar();
		exit(1);
	}

	for (int ie = 0; ie < nEns + 1; ie++)            //nEns + 1 the extra entry for ensemble mean
	{	
		for (int is = 0; is < 6; is++)
			stateVS[is][ie] = statev[is];
//#*$8.8.16 7th state for snow surface temp 
		stateVS[6][ie] = taveprevday[nstepinaDay - 1];
		//8th for snow ave. temp Tave
		stateVS[7][ie] = tsprevday[nstepinaDay - 1];
		//Rain+Melt = SWIT 9.14.17, not really a state but used in ensemble streamflow forecast
		//set to 0 here--a one time step value, 
		stateVS[8][ie] = 0.0;
		
		cumPV[ie] = cumP;
		cumEsV[ie] = cumEs;
		cumEcV[ie] = cumEc;                // Evaporation from canopy
		cumMrV[ie] = cumMr;                 // canopy melt not added
		cumGmV[ie] = cumGm;             //  Cumulative glacier melt
		cumEgV[ie] = cumEg;
		dStorageV[ie] = dStorage;
		errMBV[ie] = errMB;
	}
	//to keep track of states not updaated
	for (int is = 0; is < 6; is++)
		stateVS0[is] = statev[is];
	//8.8.16 7th state for snow surface temp 
	stateVS0[6] = taveprevday[nstepinaDay - 1]; 
	stateVS0[7] = tsprevday[nstepinaDay - 1];
	//Rain+Melt = SWIT 9.14.17, not really a state but used in ensemble streamflow forecast
	//set to 0 here--a one time step value, 
	stateVS0[8] = 0.0;

	return;
}
// function  to read forcing / weather variables control file
void uebCellDA::readInputForContr(const char* inputconFile)
{
	std::ifstream pinFile(inputconFile);
	char headerLine[256];
	//istringstream valueLine;	
	pinFile.getline(headerLine, 256);   //skip header	
	for (int i = 0; i<13; i++)
	{
		pinFile.getline(headerLine, 256, ':');
		sscanf(headerLine, "%s ", &infrContArr[i].infName);
		pinFile.getline(headerLine, 256, '\n');
		pinFile.getline(headerLine, 256, '\n');
		sscanf(headerLine, "%d ", &infrContArr[i].infType);
		//headerLine[0] = 0;
		//fscanf(pinFile,"%d\n",&svArr[i].svType);
		switch (infrContArr[i].infType)
		{
		case -1:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &infrContArr[i].infdefValue);
			break;
		case 0:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%s ", &infrContArr[i].infFile);
			break;
		case 1:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%s %s %s %d", &infrContArr[i].infFile, &infrContArr[i].infvarName, &infrContArr[i].inftimeVar, &infrContArr[i].numNcfiles);
			break;
		case 2:
			pinFile.getline(headerLine, 256, '\n');
			sscanf(headerLine, "%f ", &infrContArr[i].infdefValue);
			break;
		default:
			std::cout << "Wrong input/forcing type; has to be -1 (compute by the model), 2 (single value) , 0 (time-series text file) or 1 (3D netcdf)" << std::endl;
			std::cout << "Using default value..." << std::endl;
			break; //exit(1); 
		}
		//i++;
		//headerLine[0] = 0;
		//}
	}
	pinFile.close();
	return;
}
/*
void uebCellDA::getInpForcArr(int numNc[13], float*** RegArray[13], float &tcorVar, int ncTotaltimestep[13], MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	for (int it = 0; it < 13; it++)
	{
		if (infrContArr[it].infType == 0)
		{
			// for time series from text file read once ---- outside of this function
		}
		else if (infrContArr[it].infType == 2 || infrContArr[it].infType == -1)
		{
			// use default value or compute internally		
		}
		else if (infrContArr[it].infType == 1)     // == 0
		{
			tEnd = ncReadStart[it] + 24;          //read  24 previously 
			//offSet = 1;    // uebCellDAY*dimLen2*numTotalTimeSteps + uebCellDAX*numTotalTimeSteps;
			int retvalue = 0;			
			//read 3D netcdf (regridded array processed by uebInputs)
			char numtoStr[256];
			sprintf(numtoStr, "%d", numNc[it]);
			char tsInputfile[256];
			strcpy(tsInputfile, infrContArr[it].infFile);
			strcat(tsInputfile, numtoStr);
			strcat(tsInputfile, ".nc");
			//std::cout<<"%s\n",tsInputfile);
			//clear existing memory RegArray[it] before passing to this function // delete[] RegArray[it];
			readNC_yxSlub(tsInputfile, infrContArr[it].infvarName, infrContArr[it].inftimeVar, ncReadStart[it], tEnd, RegArray[it], ncTotaltimestep[it], numNc[it], inpComm, inpInfo);
			//retvalue = readNC_yxSlub_givenT(tsInputfile, infrContArr[it].infvarName, infrContArr[it].inftimeVar, ncReadStart[it], RegArray[it], tcorVar, numNc[it], inpComm, inpInfo);
			//startIndex[it] = 0;
			//endIndex[it] = ncTotaltimestep[it];
			//std::cout << "nc time = " << ncNtimestes[it][numNc];			
		}
	}	
}
*/
void uebCellDA::updateInpForcArr(float*** RegArray[13], int ncTotaltimestep[13])
{
	setInpForcArr(0, RegArray[0], PrecArr, ncTotaltimestep[0]);
	setInpForcArr(1, RegArray[1], TempArr, ncTotaltimestep[1]);
	setInpForcArr(2, RegArray[2], TaminArr, ncTotaltimestep[2]);
	setInpForcArr(3, RegArray[3], TamaxArr, ncTotaltimestep[3]);
	setInpForcArr(4, RegArray[4], WindspArr, ncTotaltimestep[4]);
	setInpForcArr(5, RegArray[5], RhArr, ncTotaltimestep[5]);
	setInpForcArr(6, RegArray[6], VpArr, ncTotaltimestep[6]);
	setInpForcArr(7, RegArray[7], ApresArr, ncTotaltimestep[7]);
	setInpForcArr(8, RegArray[8], SradArr, ncTotaltimestep[8]);
	setInpForcArr(9, RegArray[9], LradArr, ncTotaltimestep[9]);
	setInpForcArr(10, RegArray[10], NradArr, ncTotaltimestep[10]);
	setInpForcArr(11, RegArray[11], QgArr, ncTotaltimestep[11]);
	setInpForcArr(12, RegArray[12], SnowalbArr, ncTotaltimestep[12]);
}

void uebCellDA::setInpForcArr(int it, float ***inArray, float* forcArr, int ncTotaltimestepit)
{
	//need to call each variable array as each array has to be copied separately to device array in cuda
	int tsLength = 24; //default length	
	if (infrContArr[it].infType == 0)
	{	
		if (ncTotaltimestepit - startIndex[it] < tsLength)
			tsLength = ncTotaltimestepit - startIndex[it];  //make sure not to go out of array bound
		if (numSimTimeSteps > tsLength)
			numSimTimeSteps = tsLength;            // use the smallest number of sim time steps based on available data

		for (int i = 0; i < tsLength; i++)
			forcArr[i] = inArray[0][0][startIndex[it] + i];
		startIndex[it] += tsLength;
	}
	else if (infrContArr[it].infType == 2 || infrContArr[it].infType == -1)
	{
		// use default value or compute internally		
	}
	else if (infrContArr[it].infType == 1)     // == 0
	{
		if (ncTotaltimestepit - startIndex[it] < tsLength)
			tsLength = ncTotaltimestepit - startIndex[it];  //make sure not to go out of array bound
		if (numSimTimeSteps > tsLength)
			numSimTimeSteps = tsLength;            // use the smallest number of time steps
		
		for (int i = 0; i < tsLength; i++)
			forcArr[i] = inArray[startIndex[it] + i][uebCellDAY][uebCellDAX];
		startIndex[it] = 0;       // +tsLength;
	}

	//need to call each variable array as each array has to be copied separately to device array in cuda
	/*if (infrContArr[it].infType == 0)
	{
	forcArr = inArray[0][startIndex[it]];
	startIndex[it] = startIndex[it] + 1;
	}
	else if (infrContArr[it].infType == 2 || infrContArr[it].infType == -1)
	{
	// use default value or compute internally
	}
	else if (infrContArr[it].infType == 1)     // == 0
	forcArr = inArray[uebCellDAY][uebCellDAX];*/
		
}
// read input text file and record datetime, value pair --skip no data, get no data value from file
void  uebCellDA::readTStextFileTimeValPair(const char* inforcFile, std::vector<std::pair<double, float> > &tvar_in, int &nrecords)
{
	std::ifstream inputFile(inforcFile, std::ios::in);
	if (!inputFile)
	{
		std::cout << "Error opening file: " << inforcFile << std::endl;
		return;
	}
	nrecords = 0;
	float noDataV = -9999;
	char commentLine[256];                    //string to read header line
	inputFile.getline(commentLine, 256, '\n');
	sscanf(commentLine, "%f \n", &noDataV); //   get no data value from file

	inputFile.getline(commentLine, 256, '\n');  //skip header line 
	int Year, Month, Day;
	double Hour, DTimeV;
	float Value;
	while (!inputFile.eof())
	{
		commentLine[0] = ' ';
		inputFile.getline(commentLine, 256, '\n');
		if (commentLine[0] != ' ') {  //condition to make sure empty line is not read; 
			sscanf(commentLine, "%d %d %d %lf %f \n", &Year, &Month, &Day, &Hour, &Value); //
			//std::cout << " hour " << Hour;
			if (fabs(Value - noDataV) > 0.1) {             //only copy data that is not no-data
				DTimeV = julian(Year, Month, Day, Hour);
				//std::cout << " hour julian " << std::setprecision(15)<< DTimeV;
				tvar_in.push_back(std::make_pair(DTimeV, Value));
				++nrecords;
			}
		}
	}//while
	inputFile.close();

	return;
}// __host__ __device__  void    
//print all output values at a point
void uebCellDA::printPointOutputs(const char* outFileName)
{
	FILE* outFile = fopen(outFileName,"a");             //can write multiple times appending at the end
	//for (int istep = 0; istep < numSimTimeSteps - 1; istep++){         //-2 to be safe againts ceil( ) in computeModelDateTime()
	fprintf(outFile, "\n %d %d %d %8.3f ", (int)OutVarValues[0], (int)OutVarValues[1], (int)OutVarValues[2], OutVarValues[3]);
	for (int vnum = 4; vnum <71; vnum++)
		fprintf(outFile, " %16.4f ", OutVarValues[vnum]);
	//}
	fclose(outFile);
}
//print values at a point for degugging
void uebCellDA::printDebugOutputs()
{
	char testPrint[256];
	char ind[256];
	strcpy(testPrint, "ZTest");
	sprintf(ind, "%d", uebCellDAY);
	strcat(testPrint, ind);
	strcat(testPrint, "_");
	sprintf(ind, "%d", uebCellDAX);
	strcat(testPrint, ind);
	strcat(testPrint, ".txt");	
	FILE* outFile = fopen(testPrint, "a");
	//for (int istep = 0; istep < numSimTimeSteps - 1; istep++){         //-2 to be safe againts ceil( ) in computeModelDateTime()
	fprintf(outFile, "\n %d %d %d %8.3f ", (int)OutVarValues[0], (int)OutVarValues[1], (int)OutVarValues[2], OutVarValues[3]);
	for (int vnum = 4; vnum < 71; vnum++)
		fprintf(outFile, " %16.4f ", OutVarValues[vnum]);
	//}
	fclose(outFile);
}
//print SWE (snow water equivalent), Us (Energy content), P(recipitation) and Ta(Temperature) at a point
void uebCellDA::printSampleOutputs(const char* outFileName)
{
	FILE* outFile = fopen(outFileName,"a");
	//for (int istep = 0; istep < numSimTimeSteps - 1; istep++)         //-2 to be safe againts ceil( ) in computeModelDateTime()	
	fprintf(outFile, "\n %d %d %d %8.3f %16.4f %16.4f %16.4f %16.4f %16.4f %16.4f", (int)OutVarValues[0], (int)OutVarValues[1], (int)OutVarValues[2],
		OutVarValues[3], OutVarValues[12], OutVarValues[13], OutVarValues[16], OutVarValues[17], OutVarValues[18], OutVarValues[19]);
	fclose(outFile);
}
__host__ __device__ 
int uebCellDA::findMax(int a, int  b)
{
	return (a>b)?a:b;	

}
__host__ __device__ 
int uebCellDA::findMin(int a, int b)
{
	return (a<b)?a:b;
}
__host__ __device__ 
float uebCellDA::findMax(float a, float  b)
{
	return (a>b)?a:b;
}
__host__ __device__ 
float uebCellDA::findMin(float a, float b)
{
	return (a<b)?a:b;
}
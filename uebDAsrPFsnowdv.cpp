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

//copied (and modified) from snowdv.f90  
#include "uebDAsrPFuebpgdecls.h"
#include "Eigen/Dense"
using namespace Eigen;

//this uses class forcing arrays
/*__host__ __device__ void uebCellDA::runUEB(int dimLen2)
{
	//dimLen2
	runUEB();
}*/
__host__ __device__ 
int uebCellDA::getforcOffset(int ioffst, int dimLen2)
{
	if (ioffst == 1)
		return uebCellDAY*dimLen2*numTotalTimeSteps + uebCellDAX*numTotalTimeSteps;
	else
		return 0;
}
__host__ __device__ 
void uebCellDA::setForcingAndRadiationParamterization()
{
	//int indx = blockIdx.x*blockDim.x + threadIdx.x;
	// UTC to local time conversion
	calendardate(currentModelDateTime, Year, Month, Day, dHour);
	UTCHour = dHour - UTCOffset;
	OHour = UTCHour + lon / 15.0;
	UTCJulDat = julian(Year, Month, Day, OHour);
	calendardate(UTCJulDat, MYear, MMonth, MDay, MHour);
	fHour = (float)MHour;
	fModeldt = (float)modelDT;

	//copy data for observed for da ----7.19.16 till better way found
	//      Map from wrapper input variables to UEB variables				
	P = PrecArr[istep];                               // / 24000;   #12.19.14 --Daymet prcp in mm/day
	V = WindspArr[istep];
	Ta = TempArr[istep];

	if (infrContArr[2].infType == 0 || infrContArr[2].infType == 1)
	{
		Tmin = TaminArr[istep];
		Tmax = TamaxArr[istep];
	}
	else
	{
		//get min max temp
		Tmax = TempArr[istep];
		Tmin = TempArr[istep];
		//get max/min temperature during the day 
		nb = (dHour - 0) / modelDT;		 //number of time steps before current time within same day
		nf = (24 - dHour) / modelDT;      //no of time steps after current time within the same day
										  //#_TBC 9.13.13 look for better method for the following
		if (dHour > 23)                  //to take care of hour 24, <=>0hr
		{
			nb = 0;
			nf = 24 / modelDT;
		}
		nbstart = findMax(istep - nb, 0);  //to guard against going out of lower bound near start time when the start time is not 0 hr (istep < nb )
		nfend = findMin(istep + nf, numSimTimeSteps - 1); //don't go out of upper limit

		for (int it = nbstart; it < nfend; ++it)
		{
			if (TempArr[it] <= Tmin)
				Tmin = TempArr[it];
			if (TempArr[it] >= Tmax)
				Tmax = TempArr[it];
		}
	}
	Trange = Tmax - Tmin;
	if (snowdgtvariteflag == 1)
		std::cout<<"Temperature range "<<Trange<<std::endl;
	if (Trange <= 1)
	{
		if (snowdgtvariteflag == 1)
		{
			std::cout << "Input Diurnal temperature range "<<Trange<<" is less than or equal to 1 which is unrealistic " << std::endl;
			std::cout << "Diurnal temperature range is assumed as 8 degree celsius for grid cell y= " << uebCellDAY << " x= " << uebCellDAX << " on ";
			std::cout << Year << " " << Month << " " << Day << std::endl;
		}
		Trange = 8.0;
	}
	//		 Flag to control radiation (irad)
	//!     0 is no measurements - radiation estimated from diurnal temperature range
	//!     1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
	//!     2 is incoming shortwave and longwave radiation read from file (measured)
	//!     3 is net radiation read from file (measured)
	switch (irad)
	{
	case 0:
		Qsiobs = infrContArr[8].infdefValue;
		Qli = infrContArr[9].infdefValue;
		Qnetob = infrContArr[10].infdefValue;
		break;
	case 1:
		Qsiobs = SradArr[istep] ;                         // *3.6; // Daymet srad in W/m^2
		Qli = infrContArr[9].infdefValue;
		Qnetob = infrContArr[10].infdefValue;
		break;
	case 2:
		Qsiobs = SradArr[istep] ;                         // *3.6; // Daymet srad in W/m^2
		Qli = LradArr[istep] ;
		Qnetob = infrContArr[10].infdefValue;
		break;
	case 3:
		Qsiobs = infrContArr[8].infdefValue;
		Qli = infrContArr[9].infdefValue;
		Qnetob = NradArr[istep] ;
		break;
	default:
		std::cout << " The radiation flag is not the right number; must be between 0 and 3" << std::endl;
		getchar();
		break;
	}
	//atm. pressure from netcdf 10.30.13    //this needs revision 		//####TBC_6.20.13
	if (infrContArr[7].infType == 2)
		sitev[1] = infrContArr[7].infdefValue;
	else
		sitev[1] = ApresArr[istep] ;
	//this needs revision 		//####TBC_6.20.13
	if (infrContArr[11].infType == 2)
		Qg = infrContArr[11].infdefValue;
	else
		Qg = QgArr[istep] ;
	//!     Flag to control albedo (ireadalb)  				 
	if (infrContArr[12].infType == 2)
		Snowalb = infrContArr[12].infdefValue;
	else
		Snowalb = SnowalbArr[istep] ;
	//12.18.14 Vapor pressure of air
	if (infrContArr[6].infType == 2)
		Vp = infrContArr[6].infdefValue;
	else
		Vp = VpArr[istep] ;
	//relative humidity computed or read from file
	//#12.18.14 needs revision
	if (infrContArr[5].infType == 2)
	{
		RH = infrContArr[5].infdefValue;
	}
	else if (infrContArr[5].infType == -1)          //RH computed internally 
	{
		float eSat = svp(Ta);            //vapor pressure over water or ice depending on Ta       611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
		RH = Vp / eSat;
	}
	else
		RH = RhArr[istep] ;
	if (RH > 1)
	{
		if (snowdgtvariteflag == 1)
			std::cout<<"relative humidity >= 1 at time step "<<istep<<std::endl;
		RH = 0.99;
	}

	//  Below is code from point UEB 
	sitev[2] = Qg;
	Inpt[0] = Ta;
	Inpt[1] = P;
	Inpt[2] = V;
	Inpt[3] = RH;
	Inpt[6] = Qnetob;

	//Radiation Input Parameterization  
	hyri(MYear, MMonth, MDay, fHour, fModeldt, slope, azi, lat, HRI, cosZen);
	Inpt[7] = cosZen;
	if (irad <= 2)
	{
		atf(atff, Trange, Month, dtbar, bca, bcc);
		// We found that Model reanalysis and dowscaled data may produce some unreasonably negative solar radiation. this is simply bad data and it is generally better policy to try to give a model good data. 
		// If this is not possible, then the UEB checks will avoid the model doing anything too bad, it handles negative solar radiation in following way:
		// "no data in radiation would be to switch to the temperature method just for time steps when radiation is negative." 

		if (irad == 0 || Qsiobs < 0)     //  For cases where input is strictly negative we calculate QSI from HRI and air temp range.  This covers the case of missing data being flagged with negative number, i.e. -9999.                 
		{
			Inpt[4] = atff* Io *HRI;
			cloud(as, bs, atff, cf);   // For cloudiness fraction
		}
		else   // Here incoming solar is input
		{
			//      Need to call HYRI for horizontal surface to perform horizontal measurement adjustment
			hyri(MYear, MMonth, MDay, fHour, fModeldt, 0.0, azi, lat, HRI0, cosZen);
			//      If HRI0 is 0 the sun should have set so QSIOBS should be 0.  If it is
			//      not it indicates a potential measurement problem. i.e. moonshine
			if (HRI0 > 0)
			{
				//std::cout<<Qsiobs;
				atfimplied = findMin(Qsiobs / (HRI0*Io), 0.9); // To avoid unreasonably large radiation when HRI0 is small
				Inpt[4] = atfimplied * HRI * Io;
			}
			else
			{
				Inpt[4] = Qsiobs;
				if (Qsiobs != 0)
				{
					if (radwarnflag < 3)   //leave this warning only three times--enough to alert to non- -ve night time solar rad
					{
						std::cout << "Warning: you have nonzero nightime incident radiation of " << Qsiobs << std::endl;
						std::cout << "at date " << Year << "   " << Month << "   " << Day << "     " << dHour << std::endl;
						++radwarnflag;
					}
				}
			}
			cloud(as, bs, atff, cf);   // For cloudiness fraction  This is more theoretically correct
		}
		if (irad < 2)
		{
			qlif(Ta, RH, T_k, SB_c, Ema, Eacl, cf, QLif);
			Inpt[5] = QLif;
		}
		else
		{
			Ema = -9999;  //  These values are not evaluated but may need to be written out so are assigned for completeness
			Eacl = -9999;
			Inpt[5] = Qli;
		}
		iradfl = 0;
	}   // Long wave or shortwave either measured and calculated
	else
	{
		iradfl = 1;                    // This case is when given IRAD =3 (From Net Radiation)  
		Inpt[6] = Qnetob;
	}

	//      set control flags
	iflag[0] = iradfl;   // radiation [0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7]
	//  In the code above radiation inputs were either computed or read from input files
	iflag[1] = 0;        // no 0 [/yes 1] printing
	//iflag[2] = outFile;        // Output unit to which to print
	if (ireadalb == 0)
		iflag[3] = 1;        // Albedo Calculation [a value 1 means albedo is calculated, otherwise statev[3] is albedo
	else
	{
		iflag[3] = 0;
		statev[2] = Snowalb;
	}
	
	mtime[0] = Year;
	mtime[1] = Month;
	mtime[2] = Day;
	mtime[3] = dHour;

	return;
}
//11.27.17 for RDHM only one time step considered at a time
//inputs are P, Ta, Tmin(daily), Tmax(daily), VP (humidity), V
__host__ __device__ 
void uebCellDA::setForcingAndRadiationParamterization(float P, float Ta, float Tmin, float Tmax, float VP, float V)
{
	//int indx = blockIdx.x*blockDim.x + threadIdx.x;
	// UTC to local time conversion
	calendardate(currentModelDateTime, Year, Month, Day, dHour);
	UTCHour = dHour - UTCOffset;
	OHour = UTCHour + lon / 15.0;
	UTCJulDat = julian(Year, Month, Day, OHour);
	calendardate(UTCJulDat, MYear, MMonth, MDay, MHour);
	fHour = (float)MHour;
	fModeldt = (float)modelDT;

	//daily temp range
	Trange = Tmax - Tmin;
	if (snowdgtvariteflag == 1)
		std::cout << "Temperature range " << Trange << std::endl;
	if (Trange <= 1)
	{
		if (snowdgtvariteflag == 1)
		{
			std::cout << "Input Diurnal temperature range " << Trange << " is less than or equal to 1 which is unrealistic " << std::endl;
			std::cout << "Diurnal temperature range is assumed as 8 degree celsius for grid cell y= " << uebCellDAY << " x= " << uebCellDAX << " on ";
			std::cout << Year << " " << Month << " " << Day << std::endl;
		}
		Trange = 8.0;
	}

//*****11.27.17: 
	//rad from daily temperature profile
	//atm. pressure from sitevar
	// Qt from sitevars (set as QgArr[0])
	// Snowalb from sitevars (set as SnowalbArr[0])
	Qg = QgArr[0];
	sitev[2] = Qg;			 
	Snowalb = SnowalbArr[0];

	//11.27.17 relative humidity computed from Vapor pressure of air //#12.18.14 needs revision
	//RH computed internally 
	float eSat = svp(Ta);            //vapor pressure over water or ice depending on Ta          611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
	RH = VP / eSat;
	if (RH > 1)
	{
		if (snowdgtvariteflag == 1)
			std::cout << "relative humidity >= 1 at time step " << istep << std::endl;
		RH = 0.99;
	}

	//  Below is code from point UEB 
	Inpt[0] = Ta;
	Inpt[1] = P;
	Inpt[2] = V;
	Inpt[3] = RH;

	//Radiation Input Parameterization  from Ta
	hyri(MYear, MMonth, MDay, fHour, fModeldt, slope, azi, lat, HRI, cosZen);
	Inpt[7] = cosZen;
	atf(atff, Trange, Month, dtbar, bca, bcc);
	// We found that Model reanalysis and dowscaled data may produce some unreasonably negative solar radiation. this is simply bad data and it is generally better policy to try to give a model good data. 
	// If this is not possible, then the UEB checks will avoid the model doing anything too bad, it handles negative solar radiation in following way:
	// "no data in radiation would be to switch to the temperature method just for time steps when radiation is negative." 
	Inpt[4] = atff* Io *HRI;
	cloud(as, bs, atff, cf);   // For cloudiness fraction

	qlif(Ta, RH, T_k, SB_c, Ema, Eacl, cf, QLif);
	Inpt[5] = QLif;
	iradfl = 0;
	Inpt[6] = -9999;

	//      set control flags
	iflag[0] = iradfl;   // radiation [0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7]
						 //  In the code above radiation inputs were either computed or read from input files
	iflag[1] = 0;        // no 0 [/yes 1] printing
						 //iflag[2] = outFile;        // Output unit to which to print
	if (ireadalb == 0)
		iflag[3] = 1;        // Albedo Calculation [a value 1 means albedo is calculated, otherwise statev[3] is albedo
	else
	{
		iflag[3] = 0;
		statev[2] = Snowalb;
	}

	mtime[0] = Year;
	mtime[1] = Month;
	mtime[2] = Day;
	mtime[3] = dHour;

	return;
}
__host__ __device__
void uebCellDA::updateSimTime()
{
	istep++;
	//if (istep >= numSimTimeSteps) istep = 0;  //rewind istep to read data for next batch of sim

	UPDATEtime(Year, Month, Day, dHour, modelDT);
	currentModelDateTime = julian(Year, Month, Day, dHour);
	//copy next time step
	modelStartDate[0] = Year;
	modelStartDate[1] = Month;
	modelStartDate[2] = Day;
	modelStartHour = dHour;
	return;
}
//11.27.17--for RDHM single time step considered at a time 
//__host__ __device__
/*void uebCellDA::updateSimTime(int Year, int Month, int Day, double dHour)
{
istep++;
//if (istep >= numSimTimeSteps) istep = 0;  //

//2.20.18 model time dictated by the RDHM driver
//UPDATEtime(Year, Month, Day, dHour, modelDT);
currentModelDateTime = julian(Year, Month, Day, dHour);
//copy next time step
modelStartDate[0] = Year;
modelStartDate[1] = Month;
modelStartDate[2] = Day;
modelStartHour = dHour;
return;
}*/
//run ueb 'ordinarily'--no ensemble/ no data assimilation
__host__ __device__
void uebCellDA::runUEB()
{
	//
	SNOWUEB2();     
		  
	dStorage = statev[1]-Ws1+ statev[3]-Wc1;
	errMB= cumP-cumMr-cumEs-cumEc-dStorage+cumGm-cumEg; 				
				
	OutVarValues[0] = Year;
	OutVarValues[1] = Month;
	OutVarValues[2] = Day;
	OutVarValues[3] = dHour;
	OutVarValues[4] = atff;
	OutVarValues[5] = HRI;
	OutVarValues[6] = Eacl;
	OutVarValues[7] = Ema;
	OutVarValues[8] = Inpt[7]; //cosZen
	OutVarValues[9] = Inpt[0];
	OutVarValues[10] = Inpt[1];
	OutVarValues[11] = Inpt[2];
	OutVarValues[12] = Inpt[3];
	OutVarValues[13] = Inpt[4];
	OutVarValues[14] = Inpt[5];
	OutVarValues[15] = Inpt[6];	
							
	for (int i=16;i<69;i++)
	{		   
		OutVarValues[i]  = OutArr[i-16];					
	}
	OutVarValues[69]  = errMB;
	OutVarValues[70] = Trange;

	if (snowdgt_outflag == 1 )        //if debug mode 
	{
		printf(" time step: %d\n", istep);
		for (int uit = 0; uit<3; uit++)
			printf(" %d   ", (int) OutVarValues[uit] );
		for(int uit = 3; uit< 71; uit++)
			printf(" %16.4f  ", OutVarValues[uit] );
		printf(" \n");				
		printf("ErrMB = %16.4f \n", errMB);			
	}
	
	//to keep track of states not updated
	for (int is = 0; is < 6; is++)
		stateVS0[is] = statev[is];
	//8.8.16 7th state for snow surface temp 
	stateVS0[6] = taveprevday[nstepinaDay - 1]; 
	stateVS0[7] = tsprevday[nstepinaDay - 1];
	//SWIT 9.14.17
	stateVS0[8] = OutArr[9];	// NOTE: unit mm/hr	

	return;
}

// 3.20.18 for ens run on device side: copy historical forc
__host__ __device__
void uebCellDA::copyEnsForcingHist(float *uebForcEnsArray, int ensNum, int cellRank, int numCells)              //, float* &ensAnomaly, float &ensMean)   //, bool NormalDist)
{   
	//float P, float Ta, float Tmin, float Tmax, float VP, float V)
	if (ensNum > 59)
	{
		std::cout << " error! ensemble size exceeds 59, which is the max allowed" << std::endl; 
		std::getchar();
		exit(1);
	}
	PrecArr[ensNum] = uebForcEnsArray[cellRank];
	TempArr[ensNum] = uebForcEnsArray[numCells + cellRank];
	TaminArr[ensNum] = uebForcEnsArray[2 * numCells + cellRank];
	TamaxArr[ensNum] = uebForcEnsArray[3 * numCells + cellRank];
	VpArr[ensNum] = uebForcEnsArray[4 * numCells + cellRank];
	WindspArr[ensNum] = uebForcEnsArray[5 * numCells + cellRank];

	return;
}
//when historical forcing--including rad parametrization
__host__ __device__
void uebCellDA::runUEBEnsembles_HistForcing(int nEns)
{
	//float P, float Ta, float Tmin, float Tmax, float VP, float V
	calendardate(currentModelDateTime, Year, Month, Day, dHour);
	UTCHour = dHour - UTCOffset;
	OHour = UTCHour + lon / 15.0;
	UTCJulDat = julian(Year, Month, Day, OHour);
	calendardate(UTCJulDat, MYear, MMonth, MDay, MHour);
	fHour = (float)MHour;
	fModeldt = (float)modelDT;

	//*****11.27.17: 
	//rad from daily temperature profile
	//atm. pressure from sitevar
	// Qt from sitevars (set as QgArr[0])
	// Snowalb from sitevars (set as SnowalbArr[0])
	Qg = QgArr[0];
	sitev[2] = Qg;
	Snowalb = SnowalbArr[0];

	//      set control flags
	iflag[0] = iradfl;   // radiation [0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7]
						 //  In the code above radiation inputs were either computed or read from input files
	iflag[1] = 0;        // no 0 [/yes 1] printing
						 //iflag[2] = outFile;        // Output unit to which to print
	if (ireadalb == 0)
		iflag[3] = 1;        // Albedo Calculation [a value 1 means albedo is calculated, otherwise statev[3] is albedo
	else
	{
		iflag[3] = 0;
		statev[2] = Snowalb;
	}

	mtime[0] = Year;
	mtime[1] = Month;
	mtime[2] = Day;
	mtime[3] = dHour;
	//outputs
	OutVarValues[0] = Year;
	OutVarValues[1] = Month;
	OutVarValues[2] = Day;
	OutVarValues[3] = dHour;
	OutVarValues[4] = atff;
	OutVarValues[5] = HRI;
	OutVarValues[6] = Eacl;
	OutVarValues[7] = Ema;
	OutVarValues[8] = Inpt[7]; //cosZen
							   //9.11.17: save ensemble mean as 'base-line' when no-da
	for (int i = 9; i < 70; i++)
	{
		OutVarValues[i] = 0.0;
	}
	OutVarValues[70] = Trange;


	for (int ie = 0; ie < nEns; ie++)
	{
		//copy ens values
		P = PrecArr[ie];
		Ta = TempArr[ie];
		Tmin = TaminArr[ie];
		Tmax = TamaxArr[ie];
		Vp = VpArr[ie];
		V = WindspArr[ie];
		
		//daily temp range		
		Trange = Tmax - Tmin;

		if (snowdgtvariteflag == 1)
			std::cout << "Temperature range " << Trange << std::endl;
		if (Trange <= 1)
		{
			if (snowdgtvariteflag == 1)
			{
				std::cout << "Input Diurnal temperature range " << Trange << " is less than or equal to 1 which is unrealistic " << std::endl;
				std::cout << "Diurnal temperature range is assumed as 8 degree celsius for grid cell y= " << uebCellDAY << " x= " << uebCellDAX << " on ";
				std::cout << Year << " " << Month << " " << Day << std::endl;
			}
			Trange = 8.0;
		}
		
		//11.27.17 relative humidity computed from Vapor pressure of air //#12.18.14 needs revision
		//RH computed internally 
		float eSat = svp(Ta);            //vapor pressure over water or ice depending on Ta          611 * exp(17.27*Ta / (Ta + 237.3)); //Pa
		
		RH = Vp / eSat;
		if (RH > 1)
		{
			if (snowdgtvariteflag == 1)
				std::cout << "relative humidity >= 1 at time step " << istep << std::endl;
			RH = 0.99;
		}

		//  Below is code from point UEB 
		Inpt[0] = Ta;
		Inpt[1] = P;
		Inpt[2] = V;
		Inpt[3] = RH;

		//Radiation Input Parameterization  from Ta
		hyri(MYear, MMonth, MDay, fHour, fModeldt, slope, azi, lat, HRI, cosZen);
		Inpt[7] = cosZen;
		atf(atff, Trange, Month, dtbar, bca, bcc);
		// We found that Model reanalysis and dowscaled data may produce some unreasonably negative solar radiation. this is simply bad data and it is generally better policy to try to give a model good data. 
		// If this is not possible, then the UEB checks will avoid the model doing anything too bad, it handles negative solar radiation in following way:
		// "no data in radiation would be to switch to the temperature method just for time steps when radiation is negative." 
		Inpt[4] = atff* Io *HRI;
		cloud(as, bs, atff, cf);   // For cloudiness fraction

		qlif(Ta, RH, T_k, SB_c, Ema, Eacl, cf, QLif);
		Inpt[5] = QLif;
		iradfl = 0;
		Inpt[6] = -9999;

		//std::cout<< "state variable from previous time step: " << std::endl;

	 //initialize states to corresponding ensemble state member from previous run
		for (int is = 0; is < 6; is++)
			statev[is] = stateVS[is][ie];
		//#*$8.8.16 7th state for snow surface temp 
		taveprevday[nstepinaDay - 1] = stateVS[6][ie];
		tsprevday[nstepinaDay - 1] = stateVS[7][ie];

		cumP = cumPV[ie];
		cumEs = cumEsV[ie];
		cumEc = cumEcV[ie];                // Evaporation from canopy
		cumMr = cumMrV[ie];               // canopy melt not added
		cumGm = cumGmV[ie];             //  Cumulative glacier melt
		cumEg = cumEgV[ie];

		//run model
		SNOWUEB2();

		// save ensemble states for next time step: some of these states are getting updated, some are not
		for (int is = 0; is < 6; is++)
			stateVS[is][ie] = statev[is];
		//8.8.16 8th state for snow surface temp
		stateVS[6][ie] = taveprevday[nstepinaDay - 1];
		stateVS[7][ie] = tsprevday[nstepinaDay - 1];
		//9.14.17 rmelt (SWIT)
		stateVS[8][ie] = OutArr[9];        //Note: unit in mm/hr

										   // accumulate for mass balance
		cumPV[ie] = cumP;
		cumEsV[ie] = cumEs;
		cumEcV[ie] = cumEc;                // Evaporation from canopy
		cumMrV[ie] = cumMr;                 // canopy melt not added
		cumGmV[ie] = cumGm;             //  Cumulative glacier melt
		cumEgV[ie] = cumEg;

		dStorageV[ie] = statev[1] - Ws1 + statev[3] - Wc1;
		errMBV[ie] = cumPV[ie] - cumMrV[ie] - cumEsV[ie] - cumEcV[ie] - dStorageV[ie] + cumGmV[ie] - cumEgV[ie];

		//save ensemble member state for update 
		// no need--take from StateVS

		if (snowdgt_outflag == 1)        //if debug mode 
		{
			printf(" time step: %d\n", istep);
			printf(" %d %d %d %8.4f   ", Year, Month, Day, dHour);
			printf(" %16.4f %16.4f %16.4f %16.4f", atff, HRI, Eacl, Ema);
			for (int uit = 0; uit < 8; uit++)
				printf(" %16.4f  ", Inpt[uit]);
			for (int uit = 0; uit < 53; uit++)
				printf(" %16.4f  ", OutArr[uit]);
			printf(" \n");
			printf("ErrMB = %16.4f \n", errMBV[ie]);
		}

		OutVarValues[9] += Inpt[0];
		OutVarValues[10] += Inpt[1];
		OutVarValues[11] += Inpt[2];
		OutVarValues[12] += Inpt[3];
		OutVarValues[13] += Inpt[4];
		OutVarValues[14] += Inpt[5];
		OutVarValues[15] += Inpt[6];
		for (int i = 16; i < 69; i++)
		{
			OutVarValues[i] += OutArr[i - 16];
		}
		OutVarValues[69] += errMBV[ie];   // errMB;		
	}
	//9.11.17 ens mean as nominal value
	for (int i = 9; i < 70; i++)
	{
		OutVarValues[i] /= nEns;
	}

	//compute and save ensemble mean
	for (int id = 0; id < 9; id++)
	{
		stateVS[id][nEns] = 0.0;
		for (int ie = 0; ie < nEns; ie++) {
			stateVS[id][nEns] += stateVS[id][ie];
		}
		stateVS[id][nEns] /= nEns;
	}
	//some of these states are getting updated, some are not
	for (int is = 0; is < 6; is++)
		statev[is] = stateVS[is][nEns];
	//8.8.16 8th state for snow surface temp
	taveprevday[nstepinaDay - 1] = stateVS[6][nEns];
	tsprevday[nstepinaDay - 1] = stateVS[7][nEns];
	OutVarValues[25] = stateVS[8][nEns];           // mm/ hr

	return;

}
// 3.20.18 for ens run on device side: copy spatially correlated random samples
__host__ __device__
void uebCellDA::copyEnsForcingMultiplier(int nEns, std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > ensForcingMultiplier, int cellRank)              //, float* &ensAnomaly, float &ensMean)   //, bool NormalDist)
{
	for (int ie = 0; ie < nEns; ie++)
	{
		this->ensForcingMultiplier[0][ie] = ensForcingMultiplier[0](cellRank, ie);    //Temperature 
		for (int iforc = 1; iforc < 6; iforc++)
			this->ensForcingMultiplier[iforc][ie] = ensForcingMultiplier[iforc](cellRank, ie);    //Inpt[1] 
	}

	return;
}
// 8.27.18 --corr among forc--for ens run on device side: copy spatially correlated random samples
__host__ __device__
void uebCellDA::copyEnsForcingMultiplier(int cellRank, int nEns, int totalGridLength, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>  ensForcingMultiplier,
	std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > ensForcingMultiplierVRH)              //, float* &ensAnomaly, float &ensMean)   //, bool NormalDist)
{
	if (istep == 0)
	{
		for (int ie = 0; ie < nEns; ie++)
		{
			this->ensForcingMultiplier[0][ie] = ensForcingMultiplierVRH[0](cellRank, ie);    //Temperature 
			this->ensForcingMultiplier[1][ie] = ensForcingMultiplier(0 * totalGridLength + cellRank, ie);    //Prec
			this->ensForcingMultiplier[2][ie] = ensForcingMultiplierVRH[1](cellRank, ie);    //windS 
			this->ensForcingMultiplier[3][ie] = ensForcingMultiplierVRH[2](cellRank, ie);    //RH 
			this->ensForcingMultiplier[4][ie] = ensForcingMultiplier(1 * totalGridLength + cellRank, ie);    //Qsi
			this->ensForcingMultiplier[5][ie] = ensForcingMultiplier(2 * totalGridLength + cellRank, ie);    //Qli
			//for (int iforc = 1; iforc < 6; iforc++)
			//	this->ensForcingMultiplier[iforc][ie] = ensForcingMultiplier(iforc * totalGridLength + cellRank, ie);    //Inpt[1] 
		}
	}
	else
	{
		for (int ie = 0; ie < nEns; ie++)
		{
			this->ensForcingMultiplier[0][ie] = corrFact1 * this->ensForcingMultiplier[0][ie] + corrFact2 * ensForcingMultiplierVRH[0](cellRank, ie);    //Temperature 
			this->ensForcingMultiplier[1][ie] = corrFact1 * this->ensForcingMultiplier[1][ie] + corrFact2 * ensForcingMultiplier(0 * totalGridLength + cellRank, ie);    //Prec
			this->ensForcingMultiplier[2][ie] = corrFact1 * this->ensForcingMultiplier[2][ie] + corrFact2 * ensForcingMultiplierVRH[1](cellRank, ie);    //windS 
			this->ensForcingMultiplier[3][ie] = corrFact1 * this->ensForcingMultiplier[3][ie] + corrFact2 * ensForcingMultiplierVRH[2](cellRank, ie);    //RH 
			this->ensForcingMultiplier[4][ie] = corrFact1 * this->ensForcingMultiplier[4][ie] + corrFact2 * ensForcingMultiplier(1 * totalGridLength + cellRank, ie);    //Qsi
			this->ensForcingMultiplier[5][ie] = corrFact1 * this->ensForcingMultiplier[5][ie] + corrFact2 * ensForcingMultiplier(2 * totalGridLength + cellRank, ie);    //Qli
			//for (int iforc = 1; iforc < 6; iforc++)
			//	this->ensForcingMultiplier[iforc][ie] = ensForcingMultiplier(iforc * totalGridLength + cellRank, ie);    //Inpt[1] 
		}
	}
	//if(istep == 2)	std::cout << " corFact = " << corrFact1 << " corFact2 = " << corrFact2 << std::endl;

	return;
}  //)
__host__ __device__
void uebCellDA::copyEnsForcingMultiplier(int cellRank, int nEns, int totalGridLength, float *dev_multivarNormalDistSamplesForc,
	float *dev_multivarNormalDistSamplesForcTVRH0, float *dev_multivarNormalDistSamplesForcTVRH1, float *dev_multivarNormalDistSamplesForcTVRH2)              //, float* &ensAnomaly, float &ensMean)   //, bool NormalDist)
{
	if (istep == 0)
	{
		for (int ie = 0; ie < nEns; ie++)
		{
			this->ensForcingMultiplier[0][ie] = dev_multivarNormalDistSamplesForcTVRH0[nEns * cellRank + ie];    //Temperature 
			this->ensForcingMultiplier[1][ie] = dev_multivarNormalDistSamplesForc[nEns * (0 * totalGridLength + cellRank) + ie];    //Prec
			this->ensForcingMultiplier[2][ie] = dev_multivarNormalDistSamplesForcTVRH1[nEns * cellRank + ie];    //windS 
			this->ensForcingMultiplier[3][ie] = dev_multivarNormalDistSamplesForcTVRH2[nEns * cellRank + ie];    //RH 
			this->ensForcingMultiplier[4][ie] = dev_multivarNormalDistSamplesForc[nEns * (1 * totalGridLength + cellRank) + ie];    //Qsi
			this->ensForcingMultiplier[5][ie] = dev_multivarNormalDistSamplesForc[nEns * (2 * totalGridLength + cellRank) + ie];    //Qli
			//for (int iforc = 1; iforc < 6; iforc++)
			//	this->ensForcingMultiplier[iforc][ie] = ensForcingMultiplier(iforc * totalGridLength + cellRank, ie);    //Inpt[1] 
		}
	}
	else
	{
		for (int ie = 0; ie < nEns; ie++)
		{
			this->ensForcingMultiplier[0][ie] = corrFact1 * this->ensForcingMultiplier[0][ie] + corrFact2 * dev_multivarNormalDistSamplesForcTVRH0[nEns * cellRank + ie];    //Temperature 
			this->ensForcingMultiplier[1][ie] = corrFact1 * this->ensForcingMultiplier[1][ie] + corrFact2 * dev_multivarNormalDistSamplesForc[nEns * (0 * totalGridLength + cellRank) + ie];    //Prec
			this->ensForcingMultiplier[2][ie] = corrFact1 * this->ensForcingMultiplier[2][ie] + corrFact2 * dev_multivarNormalDistSamplesForcTVRH1[nEns * cellRank + ie];    //windS 
			this->ensForcingMultiplier[3][ie] = corrFact1 * this->ensForcingMultiplier[3][ie] + corrFact2 * dev_multivarNormalDistSamplesForcTVRH2[nEns * cellRank + ie];    //RH 
			this->ensForcingMultiplier[4][ie] = corrFact1 * this->ensForcingMultiplier[4][ie] + corrFact2 * dev_multivarNormalDistSamplesForc[nEns * (1 * totalGridLength + cellRank) + ie];    //Qsi
			this->ensForcingMultiplier[5][ie] = corrFact1 * this->ensForcingMultiplier[5][ie] + corrFact2 * dev_multivarNormalDistSamplesForc[nEns * (2 * totalGridLength + cellRank) + ie];    //Qli
			//for (int iforc = 1; iforc < 6; iforc++) //	this->ensForcingMultiplier[iforc][ie] = ensForcingMultiplier(iforc * totalGridLength + cellRank, ie);    //Inpt[1] 
		}
	}
	//if(istep == 2)	std::cout << " corFact = " << corrFact1 << " corFact2 = " << corrFact2 << std::endl;

	return;
}
//radiation paramterization and other settings done before this---this function runs based on exising inputs held by the object
__host__ __device__
void uebCellDA::runUEBEnsembles(int nEns)              //, float* &ensAnomaly, float &ensMean)   //, bool NormalDist)
{
	// save the input forcing
	float Ptemp[6];
	for (int iforc = 0; iforc < 6; iforc++)
		Ptemp[iforc] = Inpt[iforc];               // [inputforcIndex];		    // Inpt[1];  

	//outputs
	OutVarValues[0] = Year;
	OutVarValues[1] = Month;
	OutVarValues[2] = Day;
	OutVarValues[3] = dHour;
	OutVarValues[4] = atff;
	OutVarValues[5] = HRI;
	OutVarValues[6] = Eacl;
	OutVarValues[7] = Ema;
	OutVarValues[8] = Inpt[7]; //cosZen
							   //9.11.17: save ensemble mean as 'base-line' when no-da
	for (int i = 9; i < 70; i++)
	{
		OutVarValues[i] = 0.0;
	}
	OutVarValues[70] = Trange;
	//std::cout<< "state variable from previous time step: " << std::endl;
	
	for (int ie = 0; ie < nEns; ie++)
	{
		//perturb forcing	
		Inpt[0] = (tempEnStdev * ensForcingMultiplier[0][ie]) + Ptemp[0];    //Temperature 
		for (int iforc = 1; iforc < 6; iforc++)
			Inpt[iforc] = (1.0 + (forcEnStdev * ensForcingMultiplier[iforc][ie])) * Ptemp[iforc];    //Inpt[1] = Ptemp*ep;
		/*if (istep == 0)
		{
			std::ofstream qsiSamples;
			qsiSamples.open("qSi_forc_samples.txt", std::ios::app);
			qsiSamples << " " << Inpt[4];
			qsiSamples.close();
		}*/

		//initialize states to corresponding ensemble state member from previous run
		for (int is = 0; is < 6; is++)
			statev[is] = stateVS[is][ie];
		//#*$8.8.16 7th state for snow surface temp 
		taveprevday[nstepinaDay - 1] = stateVS[6][ie];
		tsprevday[nstepinaDay - 1] = stateVS[7][ie];

		cumP = cumPV[ie];
		cumEs = cumEsV[ie];
		cumEc = cumEcV[ie];                // Evaporation from canopy
		cumMr = cumMrV[ie];               // canopy melt not added
		cumGm = cumGmV[ie];             //  Cumulative glacier melt
		cumEg = cumEgV[ie];

		//run model
		SNOWUEB2();

		// save ensemble states for next time step: some of these states are getting updated, some are not
		for (int is = 0; is < 6; is++)
			stateVS[is][ie] = statev[is];
		//8.8.16 8th state for snow surface temp
		stateVS[6][ie] = taveprevday[nstepinaDay - 1];
		stateVS[7][ie] = tsprevday[nstepinaDay - 1];
		//9.14.17 rmelt (SWIT)
		stateVS[8][ie] = OutArr[9];        //Note: unit in mm/hr

		// accumulate for mass balance
		cumPV[ie] = cumP;
		cumEsV[ie] = cumEs;
		cumEcV[ie] = cumEc;                // Evaporation from canopy
		cumMrV[ie] = cumMr;                 // canopy melt not added
		cumGmV[ie] = cumGm;             //  Cumulative glacier melt
		cumEgV[ie] = cumEg;

		dStorageV[ie] = statev[1] - Ws1 + statev[3] - Wc1;
		errMBV[ie] = cumPV[ie] - cumMrV[ie] - cumEsV[ie] - cumEcV[ie] - dStorageV[ie] + cumGmV[ie] - cumEgV[ie];

		//save ensemble member state for update 
		// no need--take from StateVS

		if (snowdgt_outflag == 1)        //if debug mode 
		{
			printf(" time step: %d\n", istep);
			printf(" %d %d %d %8.4f   ", Year, Month, Day, dHour);
			printf(" %16.4f %16.4f %16.4f %16.4f", atff, HRI, Eacl, Ema);
			for (int uit = 0; uit< 8; uit++)
				printf(" %16.4f  ", Inpt[uit]);
			for (int uit = 0; uit< 53; uit++)
				printf(" %16.4f  ", OutArr[uit]);
			printf(" \n");
			printf("ErrMB = %16.4f \n", errMBV[ie]);
		}

		OutVarValues[9] += Inpt[0];
		OutVarValues[10] += Inpt[1];
		OutVarValues[11] += Inpt[2];
		OutVarValues[12] += Inpt[3];
		OutVarValues[13] += Inpt[4];
		OutVarValues[14] += Inpt[5];
		OutVarValues[15] += Inpt[6];
		for (int i = 16; i < 69; i++)
		{
			OutVarValues[i] += OutArr[i - 16];
		}
		OutVarValues[69] += errMBV[ie];   // errMB;		
	}
	//9.11.17 ens mean as nominal value
	for (int i = 9; i < 70; i++)
	{
		OutVarValues[i] /= nEns;
	}
	//the original forc must be maintained
	for (int iforc = 0; iforc < 6; iforc++)
		Inpt[iforc] = Ptemp[iforc];

	//compute and save ensemble mean
	for (int id = 0; id < 9; id++)
	{
		stateVS[id][nEns] = 0.0;
		for (int ie = 0; ie < nEns; ie++) {
			stateVS[id][nEns] += stateVS[id][ie];
		}
		stateVS[id][nEns] /= nEns;
	}
	//some of these states are getting updated, some are not
	for (int is = 0; is < 6; is++)
		statev[is] = stateVS[is][nEns];
	//8.8.16 8th state for snow surface temp
	taveprevday[nstepinaDay - 1] = stateVS[6][nEns];
	tsprevday[nstepinaDay - 1] = stateVS[7][nEns];
	OutVarValues[25] = stateVS[8][nEns];           // mm/ hr

	return;
}
//for simulation involving da but not at current time step
__host__ __device__
void uebCellDA::setNextStepStates(int nEns)
{
	//this ensures consistency when the next time step has da ==> ensembles
	for (int is = 0; is < 9; is++)
		for (int ie = 0; ie < nEns + 1; ie++)
			stateVS[is][ie] = stateVS0[is];
	/*for (int ie = 0; ie < nEns + 1; ie++)
	{
		for (int is = 0; is < 6; is++)
			stateVS[is][ie] = statev[is];
		//8.8.16 7th state for snow surface temp 
		stateVS[6][ie] = taveprevday[nstepinaDay - 1];
		stateVS[7][ie] = tsprevday[nstepinaDay - 1];
		stateVS[8][ie] = OutArr[9];		//SWIT 9.14.17
	}*/
	return;
}
//set background states for next time step to filter updated
__host__ __device__
void uebCellDA::updateBackgroundStates(int nEns, int outStateIndex, Eigen::RowVectorXf updateStateArr)
{
	//std::cout <<std::endl<< "update state variable: " <<uebStates[is]<< std::endl;		
	for (int ie = 0; ie < nEns + 1; ie++) 
	{	
		stateVS[outStateIndex][ie] = updateStateArr(ie);
	}

	//this ensures consistency when the next time step has no da--no ensemble
	//also update the states and outputs
	switch (outStateIndex)
	{
		case 0:
			statev[outStateIndex] = updateStateArr(nEns);   //Us
			OutVarValues[16] = updateStateArr(nEns);
			break;
		case 1:
			statev[outStateIndex] = updateStateArr(nEns);   //Ws
			OutVarValues[17] = updateStateArr(nEns) * 1000.0;   // unit mm
			break;
		case 2:
			statev[outStateIndex] = updateStateArr(nEns);   //TauSn
			OutVarValues[18] = updateStateArr(nEns);
			break;
		case 3:       
			statev[outStateIndex] = updateStateArr(nEns);   //Wc
			OutVarValues[56] = updateStateArr(nEns) * 1000.0;   // unit mm
			break;
		case 4:
			statev[outStateIndex] = updateStateArr(nEns);   //refD
			OutVarValues[36] = updateStateArr(nEns) * 1000.0;   // unit mm
			break;
		case 5:
			statev[outStateIndex] = updateStateArr(nEns);   //TrefD
			OutVarValues[37] = updateStateArr(nEns) * 1000.0;   // unit mm
			break;
		case 6:
			taveprevday[nstepinaDay - 1] = updateStateArr(nEns);   //Tave
			OutVarValues[29] = updateStateArr(nEns);
			break;
		case 7:
			tsprevday[nstepinaDay - 1] = updateStateArr(nEns);  //Tsurfs
			OutVarValues[30] = updateStateArr(nEns);
			break;
		case 8:       //SWIT 9.14.17
			OutArr[9] = updateStateArr(nEns);   // unit mm/hr
			OutVarValues[25] = updateStateArr(nEns);   // unit mm/hr
			break;
		default:
		{
			std::cout << std::endl << " Error the index must be between 0 and 8" << std::endl;
			std::getchar();
		}
	}

	return;
}
//4.11.18: update by weight indices
__host__ __device__
void uebCellDA::updateBackgroundStates(int nEns, int outStateIndex, Eigen::RowVectorXf updateStateArr, std::vector<int> weightIndices)
{
	//std::cout <<std::endl<< "update state variable: " <<uebStates[is]<< std::endl;
	float stateAve = 0.0;
	for (int ie = 0; ie < nEns + 1; ie++)
	{
		stateVS[outStateIndex][ie] = updateStateArr(weightIndices[ie]);
		stateAve += updateStateArr(weightIndices[ie]);
	}
	stateAve = stateAve / (nEns + 1);

	//this ensures consistency when the next time step has no da--no ensemble
	//also update the states and outputs
	switch (outStateIndex)
	{
		case 0:
			statev[outStateIndex] = stateAve;   //Us
			OutVarValues[16] = stateAve;
			break;
		case 1:
			statev[outStateIndex] = stateAve;   //Ws
			OutVarValues[17] = stateAve * 1000.0;   // unit mm
			break;
		case 2:
			statev[outStateIndex] = stateAve;   //TauSn
			OutVarValues[18] = stateAve;
			break;
		case 3:
			statev[outStateIndex] = stateAve;   //Wc
			OutVarValues[56] = stateAve * 1000.0;   // unit mm
			break;
		case 4:
			statev[outStateIndex] = stateAve;   //refD
			OutVarValues[36] = stateAve * 1000.0;   // unit mm
			break;
		case 5:
			statev[outStateIndex] = stateAve;   //TrefD
			OutVarValues[37] = stateAve * 1000.0;   // unit mm
			break;
		case 6:
			taveprevday[nstepinaDay - 1] = stateAve;   //Tave
			OutVarValues[29] = stateAve;
			break;
		case 7:
			tsprevday[nstepinaDay - 1] = stateAve;  //Tsurfs
			OutVarValues[30] = stateAve;
			break;
		case 8:       //SWIT 9.14.17
			OutArr[9] = stateAve;   // unit mm/hr
			OutVarValues[25] = stateAve;   // unit mm/hr
			break;
		default:
			{
				std::cout << std::endl << " Error the index must be between 0 and 8" << std::endl;
				std::getchar();
			}
	}

	return;
}
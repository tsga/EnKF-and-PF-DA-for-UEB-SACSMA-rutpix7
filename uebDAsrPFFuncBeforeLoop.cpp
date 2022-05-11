/*
 * uebFuncBefore.cpp
 * called before the time loop for ueb run
 */
#include <vector>
#include <iostream>
#include <sstream>
#include "boost/lexical_cast.hpp"
#include "boost/shared_array.hpp"
#include "getSelectedVerticesInConnOrder.hpp"
#include "hlrms_interfaces.hpp"
#include "getFromDeck.hpp"
#include "makeModelDataDescription.h"

#include <boost/algorithm/string.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/format.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/fstream.hpp"
#include "getInputDeck.h"
#include "selectProperties.h"
#include "GridData1D.h"
#include "boost/filesystem/path.hpp"

 //#include "getIndBasinSizes.hpp"
#include "boost/date_time/posix_time/ptime.hpp"
#include "print1DInfo.h"

//#include "ueb.h"
#include "uebDAsrPFuebpgdecls.h"
#include "uebDAsrPFdafunctions.h"
//#include <time.h>
//#include <queue>
#pragma warning(disable : 4996)

extern "C"{
#include "models.h"
}//extern "C"

using namespace Ohd_Hydro;
using namespace std;
using namespace boost;
using namespace boost::gregorian;
using namespace boost::posix_time;

extern const Ohd_Hydro::DECK deck; 

int sacFuncBefore(PixelGraph& g);
int rutpix7FuncBefore(PixelGraph& g);

int calsacFuncBefore(PixelGraph& g);
int calrutpix7FuncBefore(PixelGraph& g);

int uebDAsrPFFuncBeforeLoop( PixelGraph& g )
{
//	cerr << "entering uebFuncBefore" << endl;

	const char* ueb_parNames[] = {
		
		"ueb_irad", //Radiation control flag (0=from ta, 1= input qsi, 2= input qsi,qli 3= input qnet)

		"ueb_ireadalb", //Albedo reading control flag (0=albedo is computed internally, 1 albedo is read)

		"ueb_tr", //Temperature above which all is rain (3 C)

		"ueb_ts", //Temperature below which all is snow (-1 C)

		"ueb_ems", //Emissivity of snow (nominally 0.99)

		"ueb_cg", //Ground heat capacity (nominally 2.09 KJ/kg/C)

		"ueb_z", //Nominal meas. heights for air temp. and humidity (2m)

		"ueb_zo", // Surface aerodynamic roughness(m)

		"ueb_rho", // Snow Density (Nominally 450 kg/m^3)

		"ueb_rhog", // Soil Density (nominally 1700 kg/m^3)

		"ueb_lc", //Liquid holding capacity of snow (0.05)

		"ueb_ks", //Snow Saturated hydraulic conductivity (20 m/hr)

		"ueb_de", //Thermally active depth of soil (0.1 m)

		"ueb_avo", //Visual new snow albedo (0.95)

		"ueb_aniro", //NIR new snow albedo (0.65)

		"ueb_lans", // The thermal conductivity of fresh (dry) snow (1.0 KJ/hr-m-K) 

		"ueb_lang", // The thermal conductivity of soil (4.0 KJ/hr-m-K) 

		"ueb_wlf", //Low frequency fluctuation in deep snow/soil layer (0.0654)

		"ueb_rd1", //Amplitude correction coefficient of heat conduction (1)

		"ueb_dnews", //The threshold depth of for new snow (0.001 m)

		"ueb_emc", //Emissivity of canopy (0.98)

		"ueb_alpha", //Scattering coefficient for solar radiation (0.5)

		"ueb_alphal", //Scattering coefficient for long wave radiation (0.0)

		"ueb_g", //leaf orientation with respect to zenith angle (0.5)

		"ueb_uc", //Unloading rate coefficient (0.004626286 Per hour) (Hedstrom and Pomeroy, 1998)

		"ueb_as", //Fraction of extraterrestrial radiation on cloudy day, Shuttleworth (1993)  

		"ueb_Bs", //(as+bs):Fraction of extraterrestrial radiation on clear day, Shuttleworth   

		"ueb_lambda", //Ratio of direct atm radiation to diffuse, worked out from Dingman   

		"ueb_rimax", //Maximum value of Richardson number for stability correction

		"ueb_wcoeff", //Wind decay coefficient for the forest

		"ueb_a", //A in Bristow-Campbell formula for atmospheric transmittance

		"ueb_c", //C in Bristow-Campbell formula for atmospheric transmittance
		
		//11.29.17 these added here to deal with -ve values
		"ueb_utcOffset",	//UTC  offset

		"ueb_ts_last", // previous (last time) step surface temperature degree celsius

		"ueb_longitude", // A 2-D grid that contains the latitude at each grid 

		"ueb_latitude", // A 2-D grid that contains the latitude at each grid point 

	};

	const char* ueb_sitestateNames[] = {

		"Us", //Energy content

		"Ws", //Snow water equivalent

		"Wc", //Snow water equivalent canopy

		"Tausn", //Snow surface dimensionless age

		"ueb_df", //Drift factor multiplier

		"ueb_apr", //Average atmospheric pressure 

		"ueb_Aep", // Albedo extinction coefficient 

		"ueb_cc", // Canopy coverage fraction   

		"ueb_hcan", // Canopy height  

		"ueb_lai", // Leaf area index

		"ueb_Sbar", // Maximum snow load held per unit branch area 

		"ueb_ycage", // Forest age flag for wind speed profile parameterization  

		"ueb_slope", // A 2-D grid that contains the slope at each grid point  

		"ueb_aspect", // A 2-D grid that contains the aspect at each grid point 

		//"ueb_latitude", // A 2-D grid that contains the latitude at each grid point   

		"ueb_subalb", // Albedo (fraction 0-1) of the substrate beneath the snow (ground, or glacier)

		"ueb_subtype", // Type of beneath snow substrate encoded as (0 = Ground/Non Glacier, 1=Clean Ice/glacier, 2= Debris covered ice/glacier, 3= Glacier snow accumulation zone)

		"ueb_gsurf", // The fraction of surface melt that runs off (e.g. from a glacier)

		"ueb_b01", // Bristow-Campbell B for January (1)

		"ueb_b02", // Bristow-Campbell B for February (2)

		"ueb_b03", // Bristow-Campbell B for March(3)

		"ueb_b04", // Bristow-Campbell B for April (4)

		"ueb_b05", // Bristow-Campbell B for may (5)

		"ueb_b06", // Bristow-Campbell B for June (6)

		"ueb_b07", // Bristow-Campbell B for July (7)

		"ueb_b08", // Bristow-Campbell B for August (8)

		"ueb_b09", // Bristow-Campbell B for September (9)

		"ueb_b10", // Bristow-Campbell B for October (10) 

		"ueb_b11", // Bristow-Campbell B for November (11)

		"ueb_b12", // Bristow-Campbell B for December (12)

		//11.29.17 moved to params to deal with -ve values
		//"ueb_ts_last", // previous (last time) step surface temperature degree celsius

		//"ueb_longitude", // A 2-D grid that contains the latitude at each grid 

	};

	const char* ueb_QgAlbNames[] = {

		"ueb_Qg",   //ground energy flux

		"ueb_SnowAlb",  // snow albedo
	};

	const char* ueb_stNames[] = {

		"Us", //Energy content
		
		"Ws", //Snow water equivalent
		
		"Wc", //Snow water equivalent canopy
	
		"Tausn", //Snow surface dimensionless age

	};

	int error;
	int npix;
	int irank;

	if (error = get_npix(g, npix)) return error;

	//if (error = save_value(g, "npix", npix)) return error;

#ifdef DEBUG_GORMS
	cerr << "     npix = " << npix << endl;
#endif //#ifdef DEBUG_GORMS
         
    std::cout << " before loop number of pixels: " << npix << std::endl;

	/*std::wcout << " calling sac for test " << std::endl;	
	error = sacFuncBefore(g);
	std::cout << error << std::endl;
	std::getchar();*/

	//start & end datetime
	int ModelStartDate[3], ModelEndDate[3]; //check this	
	int ModelStartHour, ModelEndHour;
	double ModelDtMinutes, ModelDt;
	if (error = get_start_year_month_day_hour(g, ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], ModelStartHour)) return error;
	if (error = get_end_year_month_day_hour(g, ModelEndDate[0], ModelEndDate[1], ModelEndDate[2], ModelEndHour)) return error;
	// time step 
	if (error = get_time_step_minutes(g, ModelDtMinutes)) return error;
	ModelDt = (double)ModelDtMinutes / 60.0;   //dT in hours

	//ueb parameters
	char **pUEBParInputs = new char*[44];
    for(int ic=0; ic<44; ic++)
		pUEBParInputs[ic] = new char[256];
	int numUEbparInp=44;
	if (error = get_user_data(g, &pUEBParInputs[0], numUEbparInp)) 
		return error;
	//std::cout << " user data: " << std::endl;
	//for (int is = 0; is < numUEbparInp; is++)
	//	std::cout << pUEBParInputs[is] << std::endl;

	float pUEBPar[32];
	char parName[256];
	for (size_t i = 0; i < 32; ++i)
	{
		sscanf(pUEBParInputs[i], "%s %f ",parName, &pUEBPar[i]);
		//std::cout << parName <<": "<< pUEBPar[i] << std::endl;
	}
	//11.29.17 these are handled separately because they can take -ve values	
	double ModelUTCOffset;
	float uebTsLast, uebLongitude, uebLatitude;  
	sscanf(pUEBParInputs[32], "%s %lf ", parName, &ModelUTCOffset);
	//std::cout << parName << ": " << ModelUTCOffset << std::endl;
	sscanf(pUEBParInputs[33], "%s %f ", parName, &uebTsLast);
	//std::cout << parName << ": " << uebTsLast << std::endl;
	sscanf(pUEBParInputs[34], "%s %f ", parName, &uebLongitude);
	//std::cout << parName << ": " << uebLongitude << std::endl;
	sscanf(pUEBParInputs[35], "%s %f ", parName, &uebLatitude);
	//std::cout << parName << ": " << uebLatitude << std::endl;

	//2.20.18 for da
	char daconFile[256], daoutFile[256];
	char pointInputFile[256]; //10.6.17 for point site variables at obsr. (SNOTEL) stations
	char uebDAsacrutpix7State[256]; //= { "SWE" };     //, "TSURFs" };// "Wc",  "tausn", "refDepth", "totalRefDepth",
	char stateInpFile[256], doDAssimilation[256], indaQFile[256], ueb_useRDHMUnits[256];
	sscanf(pUEBParInputs[36], "%s %s ", parName, daconFile);
	//std::cout << parName << ": " << daconFile << std::endl;
	sscanf(pUEBParInputs[37], "%s %s ", parName, daoutFile);
	//std::cout << parName << ": " << daoutFile << std::endl;
	sscanf(pUEBParInputs[38], "%s %s ", parName, pointInputFile);
	//std::cout << parName << ": " << pointInputFile << std::endl;
	sscanf(pUEBParInputs[39], "%s %s ", parName, uebDAsacrutpix7State);
	//std::cout << parName << ": " << uebDAsacrutpix7State << std::endl;
	sscanf(pUEBParInputs[40], "%s %s ", parName, stateInpFile);
	//std::cout << parName << ": " << stateInpFile << std::endl;
	sscanf(pUEBParInputs[41], "%s %s ", parName, indaQFile);
	//std::cout << parName << ": " << stateInpFile << std::endl;
	sscanf(pUEBParInputs[42], "%s %s ", parName, doDAssimilation);
	//std::cout << parName << ": " << stateInpFile << std::endl;
	sscanf(pUEBParInputs[43], "%s %s ", parName, ueb_useRDHMUnits);
	//std::cout << parName << ": " << stateInpFile << std::endl;
	//std::getchar();

	//free memory for pUEBParInputs
	for (int ic = 0; ic < 44; ic++) {
		delete[] pUEBParInputs[ic];
		pUEBParInputs[ic] = NULL;
	}
	delete[] pUEBParInputs;
	pUEBParInputs = NULL;

	bool useRDHMUnits = false;
	if (strcmp(ueb_useRDHMUnits, "True") == 0 || strcmp(ueb_useRDHMUnits, "true") == 0)
		useRDHMUnits = true;
	if (error = add_name_index(g, "useRDHMUnits")) return error;
	if (error = save_value(g, "useRDHMUnits", useRDHMUnits)) return error;

	float ModeldXY_HRAP = getFromDeck<float >(deck, "pixel-size-hrap");
	std::vector<size_t> yIndxArr(npix);
	std::vector<size_t> xIndxArr(npix);
	if (error = get_row_col(g, &yIndxArr[0], &xIndxArr[0], ModeldXY_HRAP)) return error;
	//vector of active cells
	std::vector<std::pair<int, int> > activeCells;
	for (int iy = 0; iy < npix; iy++)
		activeCells.push_back(std::make_pair(yIndxArr[iy], xIndxArr[iy])); 	
	//origin coordinates ll coordinates
	double ll_XHrap, ll_YHrap;
	if (error = get_llx(g, ll_XHrap)) return error;
	if (error = get_lly(g, ll_YHrap)) return error;
    //std::cout << std::endl<< "xorg = "<<ll_XHrap <<" yorg = "<< ll_YHrap << std::endl;
	std::ofstream gridCoordinates("gridCoords.txt");
	gridCoordinates << " yorg = " << ll_YHrap << " xorg = " << ll_XHrap << std::endl;
	gridCoordinates << " dxy = " << ModeldXY_HRAP << std::endl;
	for (int id = 0; id < npix; id++)
		gridCoordinates << activeCells[id].first << " " << activeCells[id].second << std::endl;
//## WARNING: look out for grid unit inconsistency--all have to be in HRAP //2.28.18
//TODO---what if the assimilation obs. ends before simulation time end? 
	uebEnKFDA objEnKF(npix, activeCells, ll_YHrap, ll_XHrap, (const char*)daconFile, (const char*)stateInpFile, (const char*) indaQFile, (const char*)uebDAsacrutpix7State);
	//objEnKF.initDAMatrices(activeCells, daYcorrArr, daXcorrArr, R_obsErrCov, std_norm_dist_Forc_Default, norm_dist_0Mean_Tempr, norm_dist_0Mean_Default);
	objEnKF.readDAOutputControl((const char*)daoutFile);       // , inpDaControl.daOutArr, inpDaControl.daEnsArr);
	std::cout << " number of obs points " << objEnKF.numObsPoints << std::endl;
	std::cout << " Ensemble size " << objEnKF.es_enseSize << std::endl;
	std::strcpy(objEnKF.pointInputFile, pointInputFile);
	
	// create ueb model gridcell instance and copy to arrays of grid cells
	uebCellDA objCell0(pUEBPar, ModelStartDate, ModelEndDate, (double)ModelStartHour, (double)ModelEndHour, ModelDt, ModelUTCOffset, 0, 1);   //double  UTCoffset, int inpDailyorSubd, int oStride);
	if (strcmp(doDAssimilation, "False") == 0 || strcmp(doDAssimilation, "false") == 0)
		objCell0.daAssimlate = false;
	objCell0.forcEnStdev = objEnKF.forcEnStdev;
	objCell0.tempEnStdev = objEnKF.tempEnStdev; 
	objCell0.tdecorLength = objEnKF.tdecorrLength;
	//std::vector<uebCellDA>uebCellDAArr(npix, objCell0); //
	uebCellDA *uebCellDAArr = new uebCellDA[npix]; 
	std::cout << " created uebCellDA arrays and initialized member vars; number of total cells: " << npix << std::endl; // uebCellDAArr.size() << std::endl;
	//ueb site variables
	std::vector<float> uebsiteStateArray(29 * npix);	
	std::vector<float> uebQgAlb(2 * npix);
	if (error = get_values(g, ueb_sitestateNames, 29, &uebsiteStateArray[0])) return error;	
	if (error = get_values(g, ueb_QgAlbNames, 2, &uebQgAlb[0])) return error;
	float SiteState[32];
	SiteState[14] = uebLatitude;   //11.29.17 these are handled separately because they can take -ve values 
	SiteState[30] = uebTsLast;   //11.29.17 these are handled separately because they can take -ve values 
	SiteState[31] = uebLongitude;   //11.29.17 these are handled separately because they can take -ve values 
	for (irank = 0; irank < npix; irank++)
	{		
		uebCellDAArr[irank] = objCell0;
		uebCellDAArr[irank].uebCellDAY = yIndxArr[irank];
		uebCellDAArr[irank].uebCellDAX = xIndxArr[irank];
		for (int is = 0; is < 14; is++)
		{
			SiteState[is] = uebsiteStateArray[is * npix + irank];
		}
		for (int is = 15; is < 30; is++)
		{
			SiteState[is] = uebsiteStateArray[(is-1) * npix + irank];
		}

		uebCellDAArr[irank].setSiteVars_and_Initconds(SiteState);
		// Qg and snowAlb , [24];
		uebCellDAArr[irank].QgArr[0] = uebQgAlb[0 * npix + irank];
		uebCellDAArr[irank].SnowalbArr[0] = uebQgAlb[1 * npix + irank];
		/*if (irank == 0) 
		{
			for (int isout = 0; isout < 32; isout++)
				std::cout << " site var " << isout<<": "<< SiteState[isout] << std::endl;
			//std::getchar();
		}*/
		//intialize ensemble states
		uebCellDAArr[irank].setInitialEnsembleStates(objEnKF.es_enseSize);                   //, (const char*) objEnKF.daContArr.forcName);
	}
	std::cout << " done setting site vars for model domain grids" << std::endl;
	//UEB outputs: 
	//to be used inside loop for outputs
	static char timeIntStr[10];
	if (error = get_time_int_str(g, timeIntStr)) return error;

	ModelDataDescript des_rmlt;
	des_rmlt.name = "rmlt"; //name
	des_rmlt.unit = "MM";        //unit
	des_rmlt.dimension = "L";         //dimension
	des_rmlt.dataTypeCode = "RAIM";
	des_rmlt.accPolicy = SUM;     //accumulation policy
	des_rmlt.maxMissing = 0;   //max number of missing
	des_rmlt.fillMissing = DONT_FILLMISSING;    //if fill missing
	des_rmlt.missingValue = 0.f;                 //missing value
	des_rmlt.firstMissingValue = 0.f;                 //missing value for the first time step
	des_rmlt.noDataValue = -1.f;                 // no data value
	des_rmlt.timeInterval = timeIntStr;     //time interval
	if (error = add_data_descript(g, des_rmlt)) return error;
	if (error = add_name_index(g, des_rmlt.name)) return error;

	ModelDataDescript des_xmrg;
	des_xmrg.name = "xmrg"; //name
	des_xmrg.unit = "MM";        //unit
	des_xmrg.dimension = "L";         //dimension
	des_xmrg.dataTypeCode = "MAPX";
	des_xmrg.accPolicy = SUM;     //accumulation policy
	des_xmrg.maxMissing = UNLIMITED_MISSING;   //max number of missing
	des_xmrg.fillMissing = DONT_FILLMISSING;    //if fill missing
	des_xmrg.missingValue = 0.f;                 //missing value
	des_xmrg.firstMissingValue = 0.f;                 //missing value for the first time step
	des_xmrg.noDataValue = -1.f;                 // no data value
	des_xmrg.timeInterval = timeIntStr;     //time interval
	if (error = add_data_descript(g, des_xmrg)) return error;
	if (error = add_name_index(g, des_xmrg.name)) return error;

	ModelDataDescript des_uebU;
	des_uebU.name = "uebUb"; //name
	des_uebU.unit = "KJ/m^2";        //unit
	des_uebU.dimension = "FL/L^2";         //dimension
	des_uebU.dataTypeCode = "uebUb";
	des_uebU.accPolicy = AVERAGE;     //accumulation policy
	des_uebU.maxMissing = UNLIMITED_MISSING;   //max number of missing
	des_uebU.fillMissing = DONT_FILLMISSING;    //if fill missing
	des_uebU.missingValue = -99999.f;                 //missing value
	des_uebU.firstMissingValue = -99999.f;                 //missing value for the first time step
	des_uebU.noDataValue = -99999.f;                 // no data value
	des_uebU.timeInterval = timeIntStr;     //time interval
	if (error = add_data_descript(g, des_uebU)) return error;
	if (error = add_name_index(g, des_uebU.name)) return error;

	ModelDataDescript des_uebW;
	des_uebW.name = "uebW"; //name
	des_uebW.unit = "MM";        //unit
	des_uebW.dimension = "L";         //dimension
	des_uebW.dataTypeCode = "uebW";
	des_uebW.accPolicy = AVERAGE;     //accumulation policy
	des_uebW.maxMissing = 0;   //max number of missing
	des_uebW.fillMissing = DONT_FILLMISSING;    //if fill missing
	des_uebW.missingValue = 0.f;                 //missing value
	des_uebW.firstMissingValue = 0.f;                 //missing value for the first time step
	des_uebW.noDataValue = -99999.f;                 // no data value
	des_uebW.timeInterval = timeIntStr;     //time interval
	if (error = add_data_descript(g, des_uebW)) return error;
	if (error = add_name_index(g, des_uebW.name)) return error;

	//11.30.17 This keeps compatibility to UEB C++/Fortran versions
	const char* uebVars[] = { "ueb_Year", "ueb_Month", "ueb_Day", "ueb_dHour", "ueb_atff", "ueb_HRI", "ueb_Eacl", "ueb_Ema", "ueb_cosZen", "ueb_Ta", "ueb_P", "ueb_V", "ueb_RH",
		"ueb_Qsi", "ueb_Qli", "ueb_Qnet", "ueb_Us", "ueb_SWE", "ueb_tausn", "ueb_Pr", "ueb_Ps", "ueb_Alb", "ueb_QHs", "ueb_QEs", "ueb_Es", "ueb_SWIT", "ueb_QMs", "ueb_Q",
		"ueb_FM", "ueb_Tave", "ueb_TSURFs", "ueb_cump", "ueb_cumes", "ueb_cumMr", "ueb_NetRads", "ueb_smelt", "ueb_refDepth", "ueb_totalRefDepth", "ueb_cf", "ueb_Taufb",
		"ueb_Taufd", "ueb_Qsib", "ueb_Qsid", "ueb_Taub", "ueb_Taud", "ueb_Qsns", "ueb_Qsnc", "ueb_Qlns", "ueb_Qlnc", "ueb_Vz", "ueb_Rkinsc", "ueb_Rkinc", "ueb_Inmax",
		"ueb_intc", "ueb_ieff", "ueb_Ur", "ueb_Wc", "ueb_Tc", "ueb_Tac", "ueb_QHc", "ueb_QEc", "ueb_Ec", "ueb_Qpc", "ueb_Qmc", "ueb_Mc", "ueb_FMc", "ueb_SWIGM",
		"ueb_SWISM", "ueb_SWIR", "ueb_errMB", "ueb_Trange" };
	ModelDataDescript des_uebOutput;
	des_uebOutput.name = "ueb_Year"; //name
	des_uebOutput.unit = "N/A";                //unit
	des_uebOutput.dimension = "N/A";         //dimension
	des_uebOutput.dataTypeCode = "UEBOutput";
	des_uebOutput.accPolicy = AVERAGE;     //accumulation policy
	des_uebOutput.maxMissing = UNLIMITED_MISSING;   //max number of missing
	des_uebOutput.fillMissing = DONT_FILLMISSING;    //if fill missing
	des_uebOutput.missingValue = -99999.0;                 //missing value
	des_uebOutput.firstMissingValue = -99999.0;                 //missing value for the first time step
	des_uebOutput.noDataValue = -99999.0;                 // no data value
	des_uebOutput.timeInterval = timeIntStr;     //time interval
	for (int iuout = 0; iuout < 71; iuout++)
	{
		des_uebOutput.name = uebVars[iuout];
		if (error = add_data_descript(g, des_uebOutput)) return error;
		if (error = add_name_index(g, des_uebOutput.name)) return error;
		//std::cout << uebVars[iuout] << std::endl;
	}
	std::cout << "Saved ueb output variable containers " << std::endl;
	// save uebcellobjs
	if (error = add_name_index(g, "uebCellDAObjArray")) return error;
	if (error = put_value(g, &uebCellDAArr[0], "uebCellDAObjArray")) return error;

	//window outputs
	bool windowRun = false;
	//check window input from deck
	DECK::const_iterator ditr;
	ditr = deck.find("window-in-hrap");
	if (ditr != deck.end())
	{//if
		cerr << " running in a window " << endl;
		windowRun = true;
	}
	if (windowRun)
	{
		for (irank = 0; irank < npix; irank++)
		{
			char pntoStr[256];
			sprintf(pntoStr, "%d", irank);
			char outObsPoint[256];
			strcpy(outObsPoint, "obsSNOTELpoint");    // one nc file stores arrays for all variables
			strcat(outObsPoint, pntoStr);
			strcat(outObsPoint, ".txt");
			FILE* pointoutFile = fopen((const char*)outObsPoint, "w");
			for (int vnum = 0; vnum < 4; vnum++)
				fprintf(pointoutFile, "%16s", uebVars[vnum]);            // header
			for (int vnum = 4; vnum < 71; vnum++)
				fprintf(pointoutFile, "%16s ", uebVars[vnum]);            // header
			fclose(pointoutFile);
		}
	}
	if (error = add_name_index(g, "windRun")) return error;
	if (error = save_value(g, "windRun", windowRun)) return error;

//call SAC and Rutpix7 before loop
	if (error = sacFuncBefore(g)) return error;  //
	//if (error = calsacFuncBefore(g)) return error;
	std::cout << " Done SAC setup before time loop " << std::endl;  //
	if (error = rutpix7FuncBefore(g)) return error;
	//if (error = calrutpix7FuncBefore(g)) return error;
	std::cout << " Done Rutpix7 setup before time loop " << std::endl;
	/*if (error = calsacFuncBefore(g)) return error;
	std::cout << " Done SAC setup before time loop " << std::endl;
	if (error = calrutpix7FuncBefore(g)) return error;
	std::cout << " Done Rutpix7 setup before time loop " << std::endl;*/

	//for SAC Ens
	typedef boost::property_map< PixelGraph, attribute_t >::type AttMapType;
	AttMapType  attMap = boost::get(attribute_t(), g);
	typedef boost::graph_property< PixelGraph, name_index_t >::type NameIndexType;
	NameIndexType& nameIndex = boost::get_property(g, name_index_t());
	NameIndexType::iterator namePos;
	typedef bt::graph_traits< PixelGraph >::vertex_descriptor Vertex;
	std::vector< Vertex > selectedPixels;
	try { getSelectedVerticesInConnOrder<PixelGraph, std::vector >(g, selectedPixels); }//try
	catch (std::exception const& error)
	{//catch
		cerr << " getSelectedVerticesInConnOrder error!" << endl;
		return ERROR;
	}//catch
	std::vector< Vertex >::iterator itr;

	int year, month, day, hour;
	int dtm;
	if (error = get_year_month_day_hour(g, year, month, day, hour)) return error;
	if (error = get_time_step_minutes(g, dtm)) return error;
	std::vector<float> ped(npix);
	bt::shared_array< float > pe;
	bt::shared_array< float > pe_adj;
	retrieve_value(g, "pPE", pe);
	retrieve_value(g, "pPEAdj", pe_adj);
	bt::shared_array< float > sac_par;
	retrieve_value(g, "pSacPar", sac_par);

	std::vector< SACGridDA > sacGridArray(npix);

	for (int ix = 0; ix < npix; ix++)
	{
		sacGridArray[ix].dtm = (float)dtm;
		sacGridArray[ix].dtday = (float)dtm / (24. * 60);  //convert time step in days 
		sacGridArray[ix].aesc = 0;
		//sacGridArray[ix].edmnd = ped[ix];
		for (size_t i = 0; i < 17; ++i)
			sacGridArray[ix].sacpar[i] = sac_par[ix * 17 + i];
	}	
	HRAP< float > loc;
	irank = 0;
	for (itr = selectedPixels.begin(); itr != selectedPixels.end(); ++itr)
	{//for itr for (irank = 0; irank < npix; irank++)	{
		namePos = nameIndex.find("Location");
		loc = boost::any_cast<HRAP<float>>(attMap[*itr][namePos->second]);
		sacGridArray[irank].hrapx = loc.x;
		sacGridArray[irank].hrapy = loc.y;
		irank++;
	}
	//get_ped(year, month, day, npix, dtm, pe.get(), pe_adj.get(), &ped[0]);
	int pix = 0;
	for (itr = selectedPixels.begin(); itr != selectedPixels.end(); ++itr)
	{//for itr
		get_ped(year, month, day, 1, dtm, pe.get() + 12 * pix, pe_adj.get() + 12 * pix, &sacGridArray[pix].edmnd);
		namePos = nameIndex.find("ped");
		if (namePos == nameIndex.end())
		{//if
			std::cerr << "Propert: ped" << " not found, Pixel =" << *itr << std::endl;
			return PROPERTY_NOT_FOUND;
		}//if
		attMap[*itr][namePos->second] = boost::any(sacGridArray[pix].edmnd);
		//sacGridArray[pix].aesc = 0;
		pix++;
	}
	
	if (uebCellDAArr[0].daAssimlate)
	{
		int retvalue = 0;
		std::vector<uebCellDA>daUebCellObjArr(npix, objCell0); //last objEnKF.numObsPoints for obs points
		if (objEnKF.obsOutsideWS) 
		{
			//10.6.17 for site variables at da points // slope;	 aspect; 	 cc; 	 hcan; 	 lai; 	 lat;  lon; 
			float** daSiteArr = create2DArray_Contiguous(7, objEnKF.numObsPoints);    //std::vector<dasitevar> daSiteArr;
			//read da point file
			char daPointInpFile[256];
			strcpy(daPointInpFile, pointInputFile);
			strcat(daPointInpFile, "0.nc");
			//std::cout << daPointInpFile << std::endl;
			retvalue = read_Point_SiteVars_NC((const char*)daPointInpFile, daSiteArr);		// , worldComm, worldInfo); 
			//std::cout << " created uebCellDA arrays for obs points " << std::endl;
			// 10.6.17 for da cells: first fill with default (point) values --done above--then update each site variable with snotel (high res.) values
			for (int irank = 0; irank < objEnKF.numObsPoints; irank++)
			{
				//uebCellDAArr[irank] = uebCellDAArr[0];  // objCell0;
				// use point index as y coordinate and keep x coordinate as 0
				daUebCellObjArr[irank].uebCellDAY = irank + npix;
				daUebCellObjArr[irank].uebCellDAX = 0;
				//site variables from snotel stations
				SiteState[12] = daSiteArr[0][irank];	//slope
				SiteState[13] = daSiteArr[1][irank];	//aspect
				SiteState[7] = daSiteArr[2][irank];		//cc
				SiteState[8] = daSiteArr[3][irank];		//hcan
				SiteState[9] = daSiteArr[4][irank];		//lai
				SiteState[14] = daSiteArr[5][irank];	//lat
				SiteState[31] = daSiteArr[6][irank];	//lon
				daUebCellObjArr[irank].setSiteVars_and_Initconds(SiteState);
				//**** tbc latter 12.21.18
				// Qg and snowAlb , [24];
				daUebCellObjArr[irank].QgArr[0] = uebCellDAArr[0].QgArr[0];
				daUebCellObjArr[irank].SnowalbArr[0] = uebCellDAArr[0].SnowalbArr[0];
				//intialize ensemble states
				daUebCellObjArr[irank].setInitialEnsembleStates(objEnKF.es_enseSize);                   //, (const char*) objEnKF.daContArr.forcName);
				//
				//11.23.18 each obs points da without its observation (Leave one approach)
				daUebCellObjArr[irank + objEnKF.numObsPoints].uebCellDAY = irank + npix + objEnKF.numObsPoints;
				daUebCellObjArr[irank + objEnKF.numObsPoints].uebCellDAX = 0;
				daUebCellObjArr[irank + objEnKF.numObsPoints].setSiteVars_and_Initconds(SiteState);
				daUebCellObjArr[irank + objEnKF.numObsPoints].QgArr[0] = uebCellDAArr[0].QgArr[0];
				daUebCellObjArr[irank + objEnKF.numObsPoints].SnowalbArr[0] = uebCellDAArr[0].SnowalbArr[0];
				//intialize ensemble states
				daUebCellObjArr[irank + objEnKF.numObsPoints].setInitialEnsembleStates(objEnKF.es_enseSize);                   //, (const char*) objEnKF.daContArr.forcName);			
			}
			delete2DArray_Contiguous(daSiteArr);
			//std::cout << " freed point obs site var arrays" << std::endl;			
		}
		else
		{
			for (int irank = 0; irank < objEnKF.numObsPoints; irank++)
			{
				daUebCellObjArr[irank] = uebCellDAArr[objEnKF.Hc_hgVector[irank]];  // objCell0;Hc_hgVector
																					// use point index as y coordinate and keep x coordinate as 0
				daUebCellObjArr[irank].uebCellDAY = irank + npix;
				daUebCellObjArr[irank].uebCellDAX = 0;
			}
			//11.23.18 each obs points da without its observation (Leave one approach)
			for (int irank = objEnKF.numObsPoints; irank < 2 * objEnKF.numObsPoints; irank++)
			{
				daUebCellObjArr[irank] = uebCellDAArr[objEnKF.Hc_hgVector[irank - objEnKF.numObsPoints]];  // objCell0;Hc_hgVector
				// use point index as y coordinate and keep x coordinate as 0
				daUebCellObjArr[irank].uebCellDAY = irank + npix;
				daUebCellObjArr[irank].uebCellDAX = 0;
			}			
		}
		std::cout << " done setting site vars for obs points" << std::endl;

			//10.10.17 outputs at the obs points
		for (irank = 0; irank < objEnKF.numObsPoints; irank++)
		{
			char pntoStr[256];
			sprintf(pntoStr, "%d", irank);
			char outObsPoint[256];
			strcpy(outObsPoint, "obsSNOTELpoint");    // one nc file stores arrays for all variables
			strcat(outObsPoint, pntoStr);
			strcat(outObsPoint, ".txt");
			FILE* pointoutFile = fopen((const char*)outObsPoint, "w");
			for (int vnum = 0; vnum < 4; vnum++)
				fprintf(pointoutFile, "%8s", uebVars[vnum]);            // header
			for (int vnum = 4; vnum < 71; vnum++)
				fprintf(pointoutFile, "%16s ", uebVars[vnum]);            // header
			fclose(pointoutFile);
		}
		//11.23.18 each obs points da without its observation (Leave one approach)
		for (irank = 0; irank < objEnKF.numObsPoints; irank++)
		{
			char pntoStr[256];
			sprintf(pntoStr, "%d", irank);
			char outObsPoint[256];
			strcpy(outObsPoint, "obsSNOTELpoint");    // one nc file stores arrays for all variables
			strcat(outObsPoint, pntoStr);
			strcat(outObsPoint, "_leaveOut.txt");
			FILE* pointoutFile = fopen((const char*)outObsPoint, "w");
			for (int vnum = 0; vnum < 4; vnum++)
				fprintf(pointoutFile, "%8s", uebVars[vnum]);            // header
			for (int vnum = 4; vnum < 71; vnum++)
				fprintf(pointoutFile, "%16s ", uebVars[vnum]);            // header
			fclose(pointoutFile);
		}

		const char* tNameout = "time";
		//time units
		char tunits[256];
		int hhMod = (int)floor(ModelStartHour);
		int mmMod = (int)(remainder(ModelStartHour, 1.0) * 60);
		sprintf(tunits, "hours since %d-%d-%d %d:%d:00 UTC", ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], hhMod, mmMod);
		const char* tUnitsout = tunits;
		const char* tlong_name = "time";
		const char* tcalendar = "standard";

		//##2.20.18 rdhm appears to do 'one-more' simulation at the end
		int outtSteps = objCell0.numTotalTimeSteps + 1;             //--7.20.16 in the future conside saving output every outstrid'th t-step
		int outtStride = 1, outyStep = 1, outxStep = 1;
		float* t_out = new float[outtSteps];
		for (int it = 0; it < outtSteps; ++it)
			t_out[it] = it*outtStride*ModelDt;      //in hours since model start time 
		float out_fillVal = -9999.0;
		//ens da netcdf output files---for all points including obs points 1.18.18	
		int particleSize = objEnKF.es_pfSize_Fact * (objEnKF.es_enseSize + 1);
		std::cout << " Particle size = " << particleSize << std::endl;
	    /*int numEnsout = objEnKF.daEnsArr.size();
	    for (int icout = 0; icout < numEnsout - 1; icout++)
		    retvalue = Create3DNC_uebOutputs(objEnKF.daEnsArr[icout].outfName, (const char*)objEnKF.daEnsArr[icout].symbol, (const char*)objEnKF.daEnsArr[icout].units, tNameout, tUnitsout, tlong_name,
			    tcalendar, "PointIndex", "ensembleMemberNumber", outtSteps, objEnKF.numObsPoints + npix, objEnKF.es_enseSize + 1, t_out, &out_fillVal);    // , worldComm, worldInfo);
	     //retvalue = Create2DNC_uebOutputs(objEnKF.daEnsArr[numEnsout - 1].outfName, (const char*)objEnKF.daEnsArr[numEnsout - 1].symbol, (const char*)objEnKF.daEnsArr[numEnsout - 1].units, tNameout, tUnitsout, tlong_name,
		 //tcalendar, "particleMemberNumber", outtSteps, particleSize, t_out, &out_fillVal);    // , worldComm, worldInfo);*/
		
		int numDAout = objEnKF.daOutArr.size();
		for (int icout = 0; icout < numDAout - 1; icout++)
			retvalue = Create3DNC_uebOutputs(objEnKF.daOutArr[icout].outfName, (const char*)objEnKF.daOutArr[icout].symbol, (const char*)objEnKF.daOutArr[icout].units, tNameout, tUnitsout, tlong_name,
				tcalendar, "PointIndex", "ensembleMemberNumber", outtSteps, npix + (2 * objEnKF.numObsPoints) , objEnKF.es_enseSize + 1, t_out, &out_fillVal);   //, worldComm, worldInfo);
		retvalue = Create2DNC_uebOutputs(objEnKF.daOutArr[numDAout - 1].outfName, (const char*)objEnKF.daOutArr[numDAout - 1].symbol, (const char*)objEnKF.daOutArr[numDAout - 1].units, tNameout, tUnitsout, tlong_name,
			tcalendar, "particleMemberNumber", outtSteps, particleSize, t_out, &out_fillVal);
		// save cov mtrx
		/*if (error = add_name_index(g, "ueb_P_stateCovBackground")) return error;
		if (error = put_value(g, &P_stateCovBackground[0](0, 0), "ueb_P_stateCovBackground")) return error;
		if (error = add_name_index(g, "ueb_P_stateCovBackground_Points")) return error;
		if (error = put_value(g, &P_stateCovBackground_Points[0](0, 0), "ueb_P_stateCovBackground_Points")) return error;
		std::cout << "Saved ueb cov. matrices" << std::endl;*/
		if (error = add_name_index(g, "dauebCellDAObjArray")) return error;
		if (error = put_value(g, &daUebCellObjArr[0], "daUebCellObjArray")) return error;
		

		/* SAC-SMA state parameters */
		const char* sac_stnames[] = {
			"uztwc",
			"uzfwc",
			"lztwc",
			"lzfsc",
			"lzfpc",
			"adimpc"
		};
		/* SAC-SMA state parameters , in real value , not in precentage*/
		const char* sac_real_stnames[] = {
			"real_uztwc",
			"real_uzfwc",
			"real_lztwc",
			"real_lzfsc",
			"real_lzfpc",
			"real_adimpc"
		};
		// ueb forcing
		const char* ueb_forcNames[] = {
			"uebPrec",
			"uebTair",
			"uebTamin",
			"uebTamax",
			"uebVp",
			"uebWindS"
		};

		std::vector<int> rx7ndxArray(npix);
		std::vector< std::vector<float> > rx7DepthArray(npix, std::vector<float>(particleSize));
		//std::vector< Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > sacrealStateArray(npix, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(6, particleSize));   // objEnKF.es_enseSize + 1));
		std::vector < Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > rx7AreaCiArray(npix, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(4, particleSize));	  // objEnKF.es_enseSize + 1));
		if (error = get_value(g, "ndx", &rx7ndxArray[0])) return error;

		//8.1.18 perturb initial sac-rutpix7 states
		std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > multivarNormalDistSamplesForc(6,
			Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(npix, particleSize));
		for (int iforc = 0; iforc < 6; iforc++)
		{
			multivarNormalDistSamplesForc[iforc] = objEnKF.norm_dist_SACRX.samples(particleSize);
		}
		std::ofstream stdnormalSamplestxtP("particles_perturb_mult.txt");
		stdnormalSamplestxtP << multivarNormalDistSamplesForc[0] << std::endl;
//7.30.18 this is inefficient--TBRL
//##4.9.18:	//particle weight based on norm dist
		std::vector< Eigen::RowVectorXf> partWeight(npix, Eigen::RowVectorXf(particleSize));
		for (int ie = 0; ie < particleSize; ie++)		// objEnKF.es_enseSize + 1; ie++)
		{
			//initialize particle weights
			partWeight[0](ie) = 1.0;

			//#### Save ens sac and rx7 states for next step
			irank = 0;
			for (itr = selectedPixels.begin(); itr != selectedPixels.end(); ++itr)
			{//for itr for (irank = 0; irank < npix; irank++)	{
				for (int indx = 0; indx < rx7ndxArray[irank]; indx++)
				{
					//std::wcout << " here 9" << std::endl;
					namePos = nameIndex.find("areac" + boost::lexical_cast< string >(indx + 1));
					if (namePos == nameIndex.end())
					{
						std::cout << " Cannopt save areac val: Property areac" << indx + 1 << " not found at grid cell " << irank << std::endl;
						//std::getchar();
						return PROPERTY_NOT_FOUND;
					}//if
					else
					{
						//std::wcout << " here 11" << std::endl;
						rx7AreaCiArray[irank](indx, ie) = boost::any_cast<float>(attMap[*itr][namePos->second]);	// *multivarNormalDistSamplesForc[indx](irank, ie); //8118 mult.by pert.fact
						//std::wcout << " here 8" << std::endl;
					}
				}//for indx
				 //std::wcout << " here 9" << std::endl;
				namePos = nameIndex.find("depth");
				if (namePos == nameIndex.end())
				{
					std::cout << " Cannot save var depth at grid cell " << irank << std::endl;
					//std::getchar();
					return PROPERTY_NOT_FOUND;
				}//if
				else
				{
					//std::wcout << " here 11" << std::endl;
					rx7DepthArray[irank][ie] = boost::any_cast<float>(attMap[*itr][namePos->second]);		// *multivarNormalDistSamplesForc[4](irank, ie); //8118 mult.by pert.fact
					//attMap[irank][nameIndex["areac" + boost::lexical_cast<string>(indx + 1)]] = rx7AreaCiArray[indx][irank];
				}
				//std::wcout << " here 10" << std::endl;
				 //if (error = put_values(g, 6, &sacrealStateArray[ie][0], sac_real_stnames)) return error;
				for (int irstn = 0; irstn < 6; irstn++)
				{
					//std::wcout << " here 9" << std::endl;
					namePos = nameIndex.find(sac_real_stnames[irstn]);
					if (namePos == nameIndex.end())
					{
						std::cout << " Cannot save " << sac_real_stnames[irstn] << " at grid cell " << irank << std::endl;
						//std::getchar();
						return PROPERTY_NOT_FOUND;
					}//if
					else
					{
						//std::wcout << " here 11" << std::endl;
						//sacrealStateArray[irank](irstn, ie) = boost::any_cast<float>(attMap[*itr][namePos->second]);
						sacGridArray[irank].sacstEns[ie][irstn] = boost::any_cast<float>(attMap[*itr][namePos->second]);   //) * fabs(multivarNormalDistSamplesForc[irstn](irank, ie)); //8118 mult.by pert.fact
						//attMap[irank][nameIndex["areac" + boost::lexical_cast<string>(indx + 1)]] = rx7AreaCiArray[indx][irank];
					}
				}//for 				
				
				irank++;
			}
		}
		if (error = add_name_index(g, "usrParticleWeights")) return error;
		if (error = put_value(g, &partWeight[0], "usrParticleWeights")) return error;
		if (error = add_name_index(g, "rx7DepthArray")) return error;
		if (error = put_value(g, &rx7DepthArray[0], "rx7DepthArray")) return error;
		if (error = add_name_index(g, "rx7AreaCiArray")) return error;
		if (error = put_value(g, &rx7AreaCiArray[0], "rx7AreaCiArray")) return error;
		//if (error = add_name_index(g, "sacrealStateArray")) return error;
		//if (error = put_value(g, &sacrealStateArray[0], "sacrealStateArray")) return error;

		//cout << " here -10 " << endl;
		// ueb forcing for ensemble run
		typedef boost::graph_property< PixelGraph, independent_selected_basin_t >::type IndSelBasinsType;
		IndSelBasinsType independentSelectedBasins = get_property(g, independent_selected_basin_t());
		IndSelBasinsType::iterator indItr;
		std::map< std::string, size_t > bSizes;

		if (independentSelectedBasins.empty())
		{
			throw std::runtime_error(" No independent selected basins!");
		}
		for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
		{//for indItr
			size_t basinSize;
			if (!bt::get_property(g, is_connected_domain_t()))
			{//if
				basinSize = num_vertices(g);
			}//if
			else
			{//else if
				error = getNestedBasinPixelNumber(g, *indItr, basinSize);
				if (error)
				{//if
					throw std::runtime_error("getNestedBasinPixelNumer Error!");
				}//if

			}//else if
			bSizes.insert(std::make_pair(*indItr, basinSize));
		}//for indItr
		 //cout << " here -9 " << endl;
		vector< string > inLoop;
		for (int ifn = 0; ifn < 6; ifn++)
		{
			inLoop.push_back(ueb_forcNames[ifn]);
		}
		//Ohd_Hydro::GridData1D::GridData1D(std::string outdir, std::vector< std::string > const& indIds, std::map< std::string, size_t > const& sibs,
		//	std::vector< std::string > const& vars, boost::posix_time::time_duration const& s, boost::posix_time::time_period prd, bool ignore1d);
		//cout << " here -8 " << endl;
		//time_period(ptime begin, ptime end)
		boost::posix_time::ptime startPeriodEns(date(objEnKF.ModelStartDateEns[0], objEnKF.ModelStartDateEns[1], objEnKF.ModelStartDateEns[2]), hours(objEnKF.ModelStartDateEns[3]));
		boost::posix_time::ptime endPeriodEns(date(objEnKF.ModelEndDateEns[0], objEnKF.ModelEndDateEns[1], objEnKF.ModelEndDateEns[2]), hours(objEnKF.ModelEndDateEns[3]));
		boost::posix_time::time_period time_period_Ens(startPeriodEns, endPeriodEns);
		//cout << " here -7 " << endl;
		std::string outdirEns(objEnKF.xmrg1dEnsDir);
		std::cout << " outdirEns = " << outdirEns << endl;
		GridData1D gd1dEns(outdirEns,		//any_cast<string> (deck.find("output-path")->second),   //std::string  outdir,
			get_property(g, independent_selected_basin_t()),
			bSizes,				//getIndBasinSizes(g),		//std::map< std::string, size_t > const& sibs,
			inLoop,
			get_property(g, graph_time_step_t()), //any_cast< time_duration >( deck[ "time-step" ] ),
			any_cast<time_period>(time_period_Ens),     //(deck.find("time-period")->second),				//
			any_cast<bool>(deck.find("ignore-1d-xmrg")->second));
		//cout << " here -6 " << endl;
		if (error = add_name_index(g, "Ens1dXmrgObj")) return error;
		if (error = save_value(g, "Ens1dXmrgObj", gd1dEns)) return error;
		std::cout << " Before loop 1d ens info: " << endl;
		print1DInfo(gd1dEns);
		//for base-line 1d xmrg input, to preserve the forcing at the time step
		GridData1D gd1dDetr(any_cast< string > (deck.find("output-path")->second),
			get_property(g, independent_selected_basin_t()),
			bSizes,
			inLoop,
			get_property(g, graph_time_step_t()),  //any_cast< time_duration >( deck[ "time-step" ] ),
			any_cast< time_period >(deck.find("time-period")->second),
			any_cast< bool >(deck.find("ignore-1d-xmrg")->second));
		//cout << " here -6 " << endl;
		if (error = add_name_index(g, "Detr1dXmrgObj")) return error;
		if (error = save_value(g, "Detr1dXmrgObj", gd1dDetr)) return error;
		//std::cout << " Before loop 1d xmrg (determ) info: " << endl;
		//print1DInfo(gd1dDetr);
		//exit(0);

	} //if (daAssimlate == 1)
	else
	{
		const char* tNameout = "time";
		//time units
		char tunits[256];
		int hhMod = (int)floor(ModelStartHour);
		int mmMod = (int)(remainder(ModelStartHour, 1.0) * 60);
		sprintf(tunits, "hours since %d-%d-%d %d:%d:00 UTC", ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], hhMod, mmMod);
		const char* tUnitsout = tunits;
		const char* tlong_name = "time";
		const char* tcalendar = "standard";

		//##2.20.18 rdhm appears to do 'one-more' simulation at the end
		int outtSteps = objCell0.numTotalTimeSteps + 1;             //--7.20.16 in the future conside saving output every outstrid'th t-step
		int outtStride = 1, outyStep = 1, outxStep = 1;
		float* t_out = new float[outtSteps];
		for (int it = 0; it < outtSteps; ++it)
			t_out[it] = it*outtStride*ModelDt;      //in hours since model start time 
		float out_fillVal = -9999.0;

		int particleSize = 1;
		int retvalue = 0;
		int numDAout = objEnKF.daOutArr.size();

		retvalue = Create2DNC_uebOutputs(objEnKF.daOutArr[numDAout - 1].outfName, (const char*)objEnKF.daOutArr[numDAout - 1].symbol, (const char*)objEnKF.daOutArr[numDAout - 1].units, tNameout, tUnitsout, tlong_name,
			tcalendar, "particleMemberNumber", outtSteps, particleSize, t_out, &out_fillVal);
	}
	std::cout << " size of UEB obj = " << sizeof(uebCellDA) / (1024.0 * 1024.0) << " MB " << std::endl;
	std::cout << " size of SAC obj = " << sizeof(SACGridDA) / (1024.0 * 1024.0) << " MB " << std::endl;
	if (error = add_name_index(g, "uebEnKFDAObj")) return error;
	//std::vector<uebEnKFDA>objEnKFArr(1, objEnKF);
	if (error = save_value(g, "uebEnKFDAObj", objEnKF)) return error;
	//std::cout << std::endl;
	if (error = add_name_index(g, "des_sacGridArray")) return error;
	if (error = put_value(g, &sacGridArray[0], "des_sacGridArray")) return error;
	//std::cout << "Saved ueb obj and enkf obj" << std::endl;
	std::cout << " Done ueb setup before time loop " << std::endl;
	
//	cerr << "leaving uebFuncBefore" << endl;

    return OK;
}//uebFuncBeforeLoop


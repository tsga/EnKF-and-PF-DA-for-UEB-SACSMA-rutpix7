/*
* uebDAsacrutpix7FuncInsideLoop.cpp
* loop through time and grid cells, call uebrun
*/
#include <vector>
#include <iostream>
#include <sstream>
//
#include <iterator>
#include <boost/property_map/property_map.hpp>

#include "boost/lexical_cast.hpp"
#include "boost/shared_array.hpp"
//
#include <boost/graph/adjacency_list.hpp>

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

//
#include "uebDAsrPFuebpgdecls.h"
#include "uebDAsrPFdafunctions.h"

//#include "observedValue.h"
//

//#include <time.h>
//#include <queue>
#pragma warning(disable : 4996)

extern "C" {
#include "models.h"
}//extern "C"

using namespace Ohd_Hydro;
using namespace std;
using namespace boost;
using namespace boost::gregorian;
using namespace boost::posix_time;

extern const Ohd_Hydro::DECK deck; 

int sacFuncInside_Ens(PixelGraph& g, int threadsPerBlock, std::vector<float> pcp, SACGridDA *sacGridArray, int gpu0);
int sacFuncInside(PixelGraph& g);
int rutpixFuncInside(PixelGraph& g);

//int calsacFuncInside(PixelGraph& g);
//int calsacFuncAfter(PixelGraph& g);
//int calrutpixFuncInside(PixelGraph& g);
//int calrutpixFuncAfter(PixelGraph& g);


int getQ(PixelGraph& g, std::vector<float>& qArr)
{
	typedef boost::property_map< PixelGraph, attribute_t >::type AttMapType;
	AttMapType  attMap = boost::get(attribute_t(), g);
	typedef boost::graph_property< PixelGraph, id_index_t >::type IdIndexType;
	IdIndexType& idIndex = boost::get_property(g, id_index_t());

	typedef boost::graph_property< PixelGraph, selected_id_t >::type SelBasinType;
	SelBasinType selBsns = boost::get_property(g, selected_id_t());

	typedef boost::graph_property< PixelGraph, name_index_t >::type NameIndexType;
	NameIndexType& nameIndex = boost::get_property(g, name_index_t());

	NameIndexType::iterator namePos;

	namePos = nameIndex.find("discharge");
	if (namePos == nameIndex.end())
	{
		std::cerr << "Property: discharge not found" << std::endl;
		return PROPERTY_NOT_FOUND;
	}

	for (size_t bi = 0; bi < selBsns.size(); bi++)
	{
		float val = boost::any_cast<float>(attMap[boost::vertex(idIndex[std::string(selBsns[bi].c_str())], g)][namePos->second]);
		qArr.push_back(val);
	}

	return OK;
}

int putQ(PixelGraph& g, std::vector<float> qArr)
{
	typedef boost::property_map< PixelGraph, attribute_t >::type AttMapType;
	AttMapType  attMap = boost::get(attribute_t(), g);
	typedef boost::graph_property< PixelGraph, id_index_t >::type IdIndexType;
	IdIndexType& idIndex = boost::get_property(g, id_index_t());

	typedef boost::graph_property< PixelGraph, selected_id_t >::type SelBasinType;
	SelBasinType selBsns = boost::get_property(g, selected_id_t());

	typedef boost::graph_property< PixelGraph, name_index_t >::type NameIndexType;
	NameIndexType& nameIndex = boost::get_property(g, name_index_t());

	NameIndexType::iterator namePos;

	namePos = nameIndex.find("discharge");
	if (namePos == nameIndex.end())
	{
		std::cerr << "Property: discharge not found" << std::endl;
		return PROPERTY_NOT_FOUND;
	}

	for (size_t bi = 0; bi < selBsns.size(); bi++)
	{
		attMap[boost::vertex(idIndex[std::string(selBsns[bi].c_str())], g)][namePos->second] = qArr[bi];
	}

	return OK;
}

int uebDAsrPFFuncInsideLoop(PixelGraph& g)
{
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

	int error;
	int npix, irank, iuout;

	if (error = get_npix(g, npix)) return error;
	//std::cout << " Inside time loop number of pixels: " << npix << std::endl;

#ifdef DEBUG_GORMS
	cerr << "     npix = " << npix << endl;
#endif //#ifdef DEBUG_GORMS

	//std::vector<uebEnKFDA>objEnKFArr(1); // (npix);
	//if (error = get_value(g, "uebEnKFDAObjArr", &objEnKFArr)) return error;
	uebEnKFDA objEnKFArr;
	if (error = retrieve_value(g, "uebEnKFDAObj", objEnKFArr)) return error; 
	// copy to arrays of grid cells
	uebCellDA *uebCellDAObjArr = new uebCellDA[npix];
	if (error = get_value(g, "uebCellDAObjArray", &uebCellDAObjArr[0])) return error;
	//gpu control	
	int threadsPerBlock = 255;
	SACGridDA *sacGridArray = new SACGridDA[npix];
	if (error = get_value(g, "des_sacGridArray", &sacGridArray[0])) return error;
	//std::cout << " retrieved uebcellarray, da_uebCellDAObjArray, obs uebobj objects for time step " << uebCellDAObjArr[0].istep + 1 << std::endl;
	//std::cout << " retrieved uebcellarray, da_uebCellDAObjArray, obs uebobj objects for time step " << uebCellDAObjArr[0].istep + 1 << std::endl;
	
	//model domain: connectivity or window
	bool windowRun;
	if (error = retrieve_value(g, "windRun", windowRun)) return error;
	//8.7.18 use oC or F
	bool useRDHMUnits;
	if (error = retrieve_value(g, "useRDHMUnits", useRDHMUnits)) return error;

	// tamin and tamax one value per day --for the daArray from netCDF
	int currSimDay = (int)uebCellDAObjArr[0].istep / uebCellDAObjArr[0].nstepinaDay;
	//if ((uebCellDAObjArr[0].istep % 120) == 0) 	std::cout << std::endl << " UEB time step: " << uebCellDAObjArr[0].istep <<  " current day: " << currSimDay << std::endl;

//TODO---what if the assimilation obs. ends before simulation time end? 
	//start & end datetime
	int CurrentDateTime[4]; //check this	
	if (error = get_year_month_day_hour(g, CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2], CurrentDateTime[3]))	return error;
	double ModelDtMinutes, ModelDt;  //dt
	if (error = get_time_step_minutes(g, ModelDtMinutes)) return error;
	ModelDt = (double)ModelDtMinutes / 60.0;   //dT in hours

	//std::cout << " Obs R Matrix" << std::endl;	std::cout << objEnKFArr.R_obsErrCov << std::endl;	std::getchar();
	// get xmrg forcing "uebPrec",	"uebTair",	"uebTamin", "uebTamax", "uebVp", "uebWindS"

	//float *uebForcArray = new float[6*npix];
	std::vector<float> uebForcArray(6 * npix);
	if (error = get_values(g, ueb_forcNames, 6, &uebForcArray[0])) return error; 
//std::cout << " Here -3: " << std::endl;

	//"modisDelVis"
	std::vector<float> modisDelVisArray(npix);
	if (error = get_value(g, "modisDelVis",&modisDelVisArray[0])) return error; 
//std::cout << " Here -2: " << std::endl;
	for (int ipx = 0; ipx < npix; ipx++)
		uebCellDAObjArr[ipx].modisAlbedoFact = (modisDelVisArray[ipx] > 0 && modisDelVisArray[ipx] < 100) ? modisDelVisArray[ipx] : 0.0;   //100% as upper limit for dAlb

/*	if ((uebCellDAObjArr[0].istep % 120) == 0) 
	{
		std::cout << std::endl << " UEB time step: " << uebCellDAObjArr[0].istep <<  " current day: " << currSimDay << std::endl;
		std::cout << " MODIS Albedo DelVis: " << std::endl;
		for (int ipx = 0; ipx < npix; ipx++)
			std::cout << "  " << uebCellDAObjArr[ipx].modisAlbedoFact;
		std::cout << std::endl;
	}
*/
	/*for (int ifn = 0; ifn < 6; ifn++)
	{
		for (int ix = 0; ix < npix; ix++)
			std::cout << " " << uebForcArray[ifn*npix + ix];
		std::cout << std::endl;
		std::getchar();
	}
	for (int ix = 0; ix < npix; ix++)
		std::cout << " " << uebForcArray[npix + ix];
	std::cout << std::endl;*/
//std::cout << " Here -1: " << std::endl;

//float P, float Ta, float Tmin, float Tmax, float VP, float V
	float precConvFactor = 60.f / (ModelDtMinutes * 1000.0);
	if (useRDHMUnits)
	{
		//unit conversion  v (mm) during dT(min) :  m/hr => v * 60/(dt * 1000.0)	
		std::transform(uebForcArray.begin(), uebForcArray.begin() + npix, uebForcArray.begin(), std::bind1st(std::multiplies<float>(), precConvFactor));
		//Tc = 5/9(Tf - 32) [](float inV) {return 5.0*(inV - 32.0) / 9.0;}
		std::transform(uebForcArray.begin() + npix, uebForcArray.begin() + 4 * npix, uebForcArray.begin() + npix, [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
		//std::transform(uebTaminArray.begin(), uebTaminArray.end(), uebTaminArray.begin(), [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
		//std::transform(uebTamaxArray.begin(), uebTamaxArray.end(), uebTamaxArray.begin(), [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
	}
	int nOpts = objEnKFArr.numObsPoints;
	int retvalue = 0;
	//float* uebForcArrayDaPoint = new float[6 * nOpts]; 	
	std::vector<float> uebForcArrayDaPoint(6 * nOpts);
	uebCellDA *daUebCellObjArr = NULL;  
	if (uebCellDAObjArr[0].daAssimlate)
	{
		daUebCellObjArr = new uebCellDA[npix];  // , objCell0);   // 10.6.17 cells (at obs points) for data assimilation
		if (error = get_value(g, "daUebCellObjArray", &daUebCellObjArr[0])) return error;
		//std::cout << " retrieved da_uebCellDAObjArray for time step " << uebCellDAObjArr[0].istep + 1 << std::endl;
        if (objEnKFArr.obsOutsideWS) 
		{
			//if ((uebCellDAObjArr[0].istep % 120) == 0) 	std::cout << std::endl << " obs OUTside WS" << std::endl;

			char tsInputfile[256];
			std::strcpy(tsInputfile, objEnKFArr.pointInputFile);    // one nc file stores arrays for all variables
			//char enstoStr[256];
			//sprintf(enstoStr, "%d", ensN);
			//strcat(tsInputfile, enstoStr);
			std::strcat(tsInputfile, "0.nc");
			retvalue = readNC_Array_atIndex(tsInputfile, "uebPrec", "time", uebCellDAObjArr[0].istep, &uebForcArrayDaPoint[0 * nOpts]);
			retvalue = readNC_Array_atIndex(tsInputfile, "uebTair", "time", uebCellDAObjArr[0].istep, &uebForcArrayDaPoint[1 * nOpts]);
			retvalue = readNC_Array_atIndex(tsInputfile, "uebVp", "time", uebCellDAObjArr[0].istep, &uebForcArrayDaPoint[4 * nOpts]);
			retvalue = readNC_Array_atIndex(tsInputfile, "uebWindS", "time", uebCellDAObjArr[0].istep, &uebForcArrayDaPoint[5 * nOpts]);
			// tamin and tamax one value per day
			retvalue = readNC_Array_atIndex(tsInputfile, "uebTamin", "timedaily", currSimDay, &uebForcArrayDaPoint[2 * nOpts]);
			retvalue = readNC_Array_atIndex(tsInputfile, "uebTamax", "timedaily", currSimDay, &uebForcArrayDaPoint[3 * nOpts]);
			//

			//retvalue = readNC_Array_atIndex(tsInputfile, "modisDelVis", "timedaily", currSimDay, &uebForcArrayDaPoint[6 * nOpts]);

			//std::cout << " Read forcing at obs points " <<  std::endl;
			if (useRDHMUnits)
			{
				//unit conversion  v (mm) during dT(min) :  m/hr => v * 60/(dt * 1000.0)	
				std::transform(uebForcArrayDaPoint.begin(), uebForcArrayDaPoint.begin() + nOpts, uebForcArrayDaPoint.begin(), std::bind1st(std::multiplies<float>(), precConvFactor));
				//Tc = 5/9(Tf - 32) [](float inV) {return 5.0*(inV - 32.0) / 9.0;}
				std::transform(uebForcArrayDaPoint.begin() + nOpts, uebForcArrayDaPoint.begin() + 4 * nOpts, uebForcArrayDaPoint.begin() + nOpts, [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
				//std::transform(uebTaminArray.begin(), uebTaminArray.end(), uebTaminArray.begin(), [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
				//std::transform(uebTamaxArray.begin(), uebTamaxArray.end(), uebTamaxArray.begin(), [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
			}
		}
		else
		{
			//if ((uebCellDAObjArr[0].istep % 120) == 0) 	std::cout << std::endl << " obs INside WS" << std::endl;
			for (int ist = 0; ist < 6; ist++)
			{
				for (int irank = 0; irank < nOpts; irank++)
					uebForcArrayDaPoint[ist * nOpts + irank] = uebForcArray[ist * npix + objEnKFArr.Hc_hgVector[irank]];
			}		
		}		
	}
//std::cout << " Here -2: " << std::endl;
	
	//output vars
//# 2.18.18 these are only for model domain grids
	int outvarindx = 17;
	std::vector<float> UbArray(npix);
	std::vector<float> WArray(npix);
	std::vector<float> rmltArray(npix);
	std::vector<float> xmrgArray(npix);
	//std::cout << "here -1" << std::endl;
	std::vector<float> uebOutputArray(71 * npix);
	//std::vector< std::vector<float> > uebOutputArray(71, std::vector<float>(npix));
//std::cout << "here -1" << std::endl;	
	
	//std::cout<<"proc "<<rank<<" before currDT compute"<<std::endl;
	int ueb_curModelTimeStep = uebCellDAObjArr[0].istep;
	if (uebCellDAObjArr[0].daAssimlate)
	{
		if (uebCellDAObjArr[0].updateDaArray)
		{
			//1.9.17 obs da data only in txt file format for now
			//std::cout << " updating DA obs array time step " << ueb_curModelTimeStep << std::endl;
			//objEnKFArr.updateDaArr(uebCellDAObjArr[0].startIndexDA);   // daderivedType, polyThreshold, dasnIndx, Z_obs, polyCoeff1, polyCoeff2);
			uebCellDAObjArr[0].startIndexDA = uebCellDAObjArr[0].startIndexDA + 1;

//TODONE---what if the assimilation obs. ends before simulation time end? 
			if (uebCellDAObjArr[0].startIndexDA > objEnKFArr.daTcorrArr.size())
				uebCellDAObjArr[0].startIndexDA = objEnKFArr.daTcorrArr.size();   //limit the index of da array to the size of obs t array

			//std::cout << " copied obs data" << std::endl;
			//			std::cout << "proc " << rank << " DA observed state array updated" << std::endl;
			uebCellDAObjArr[0].updateDaArray = false;   // wait until this data is used in assimilation 
			/*if (uebCellDAObjArr[0].istep < 10)
			{
				std::cout << " next DA time step " << uebCellDAObjArr[0].startIndexDA << std::endl;
				std::cout << " updated DA arr at : " << uebCellDAObjArr[0].startIndexDA - 1 << " time " << objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA - 1] << std::endl;
				std::cout << objEnKFArr.Z_obs << std::endl;
			}*/
		}
		if (uebCellDAObjArr[0].updateDaArrayQ)
		{
			//1.9.17 obs da data only in txt file format for now
			uebCellDAObjArr[0].startIndexDAQ = uebCellDAObjArr[0].startIndexDAQ + 1;

//TODONE---what if the assimilation obs. ends before simulation time end? 
			if (uebCellDAObjArr[0].startIndexDAQ > objEnKFArr.daTcorrArrQ.size())
				uebCellDAObjArr[0].startIndexDAQ = objEnKFArr.daTcorrArrQ.size();   //limit the index of da array to the size of obs t array

			//std::cout << " copied obs data" << std::endl;
			//			std::cout << "proc " << rank << " DA observed state array updated" << std::endl;
			uebCellDAObjArr[0].updateDaArrayQ = false;   // wait until this data is used in assimilatio
		}
	}
	//std::cout << std::setprecision(15) << " current model datetime = " << uebCellDAObjArr[0].currentModelDateTime << ";  data assimilation next time = " << objEnKFArr.daTcorrArr [uebCellDAObjArr[0].startIndexDA - 1] << "; 0.5DT = " << 0.5*ModelDt / 24;
	//std::cout << std::setprecision(15) << "; CurrModelDateTime - NextDaDateTime = " << fabs(objEnKFArr.daTcorrArr [uebCellDAObjArr[0].startIndexDA - 1] - uebCellDAObjArr[0].currentModelDateTime) << std::endl;
	//if (rank == 0) {

//std::cout << " Here -11 " << std::endl;
	int ensSize = objEnKFArr.es_enseSize;
	int particleSize = objEnKFArr.es_pfSize_Fact * (ensSize + 1);
//std::cout << " Particle size = " << particleSize << std::endl;
	//int numEnsout = objEnKFArr.daEnsArr.size();
	int numDAout = objEnKFArr.daOutArr.size();
	float *pxvEns = new float[npix * particleSize];
	//
	typedef boost::property_map< PixelGraph, attribute_t >::type AttMapType;
	AttMapType  attMap = boost::get(attribute_t(), g);
	typedef boost::graph_property< PixelGraph, name_index_t >::type NameIndexType;
	NameIndexType& nameIndex = boost::get_property(g, name_index_t());
	NameIndexType::iterator namePos;
	typedef bt::graph_traits< PixelGraph >::vertex_descriptor Vertex;
	std::vector< Vertex >::iterator itr;
	std::vector< Vertex > selectedPixels;
	//bt::shared_array< Vertex > selectedPixels;
	//error = retrieve_value(g, "selectedPixels", selectedPixels);
	try { getSelectedVerticesInConnOrder(g, selectedPixels); }//try
	catch (std::exception const& error)
	{//catch
		cerr << " getSelectedVerticesInConnOrder error!" << endl;
		return ERROR;
	}//catch

//cout << " here 0 " << endl;
	if (uebCellDAObjArr[0].daAssimlate)
	{
/* ============when using perturbed forcing----only for model domain grid cells--obs points use only 'obs' forc */
        if (uebCellDAObjArr[0].istep == 0)
		{
			cout << " Current datetime = " << CurrentDateTime[0] << " " << CurrentDateTime[1] << " " << CurrentDateTime[2] << " " << CurrentDateTime[3] << endl;
			cout << " Ens histYear start = " << objEnKFArr.ModelStartDateEns[0] << endl;
			cout << " Ens histYear end = " << objEnKFArr.ModelStartDateEns[0] + ensSize - 1 << endl;
			cout << " Forecast datetime = " << objEnKFArr.forecastDateTimeEns[0] << " " << objEnKFArr.forecastDateTimeEns[1] << " " 
				<< objEnKFArr.forecastDateTimeEns[2] << " " << objEnKFArr.forecastDateTimeEns[3] << endl;
			cout << " Forecast datetime jul = " << objEnKFArr.forecastDateTime<<endl;
		}
		if (uebCellDAObjArr[0].currentModelDateTime < objEnKFArr.forecastDateTime)
		{
			Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> multivarNormalDistSamplesForc(3 * (npix + nOpts), ensSize);
			multivarNormalDistSamplesForc = objEnKFArr.std_norm_dist_Forc_Default.samples(ensSize);
			std::ofstream stdnormalSamplestxt("ensTempForcing_err.txt");
			stdnormalSamplestxt << multivarNormalDistSamplesForc << std::endl;

			std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > multivarNormalDistSamplesForcTVRH(3,
				Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(npix + nOpts, ensSize));
			multivarNormalDistSamplesForcTVRH[0] = objEnKFArr.norm_dist_1Mean_Default.samples(ensSize);
			multivarNormalDistSamplesForcTVRH[1] = objEnKFArr.norm_dist_1Mean_Default.samples(ensSize);
			multivarNormalDistSamplesForcTVRH[2] = objEnKFArr.norm_dist_1Mean_Default.samples(ensSize);
			std::ofstream stdnormalSamplestxtV("ensTempForcing_errT.txt");
			stdnormalSamplestxtV << multivarNormalDistSamplesForcTVRH[0] << std::endl;
			//3.20.18 copy ensemble random error samples
			/*for (irank = 0; irank < npix; irank++)
				uebCellDAObjArr[irank].copyEnsForcingMultiplier(irank, ensSize, npix + nOpts, multivarNormalDistSamplesForc, multivarNormalDistSamplesForcTVRH);
			for (irank = 0; irank < nOpts; irank++)
				daUebCellObjArr[irank].copyEnsForcingMultiplier(irank + npix, ensSize, npix + nOpts, multivarNormalDistSamplesForc, multivarNormalDistSamplesForcTVRH);
			*/
			//cout << " here -1 " << endl;
			// perturbed forcing
			//threadsPerBlock = 128;
			//callUEBrunHostEnsemble_PertForcing(threadsPerBlock, npix, ensSize, uebCellDAObjArr, uebForcArray);   // int threadsPerBlock, int npix, int nEns, uebCellDA *uebCellDAObjArr, float *uebForcArray)  //,uebEnKFDA objEnKFArr      
            /*if (objEnKFArr.obsOutsideWS) 
                callUEBrunHostEnsemble(threadsPerBlock, npix, nOpts, ensSize, uebCellDAObjArr, &uebForcArray[0], multivarNormalDistSamplesForc, multivarNormalDistSamplesForcTVRH);
			else*/
			callUEBrunHostEnsemble(threadsPerBlock, npix, nOpts, ensSize, uebCellDAObjArr, daUebCellObjArr, &uebForcArray[0], &uebForcArrayDaPoint[0], multivarNormalDistSamplesForc, multivarNormalDistSamplesForcTVRH);
			//cout << " here -10 " << endl;
		}
		else
		{
			// ueb historical forcing for ensemble run
			typedef boost::graph_property< PixelGraph, independent_selected_basin_t >::type IndSelBasinsType;
			IndSelBasinsType independentSelectedBasins = get_property(g, independent_selected_basin_t());
			IndSelBasinsType::iterator indItr;
			//cout << " here -9 " << endl;
			vector< string > inLoop;
			for (int ifn = 0; ifn < 6; ifn++)
			{
				inLoop.push_back(ueb_forcNames[ifn]);
			}

			GridData1D gd1dEns, gd1dDetr;
			if (error = retrieve_value(g, "Ens1dXmrgObj", gd1dEns)) return error;
			if (error = retrieve_value(g, "Detr1dXmrgObj", gd1dDetr)) return error;
			//boost::posix_time::time_duration timeStep(get_property(allPixels, graph_time_step_t()));	
			//std::cout << " Inside loop 1d ens info: " << endl;
			//print1DInfo(gd1dEns);
			ModelDataDescript tempDesc;
			
			//float *uebForcEnsArray = new float[6 * npix];
			std::vector<float> uebForcEnsArray(6 * npix);

			//if (uebCellDAObjArr[0].currentModelDateTime == objEnKFArr.forecastDateTime)   //(uebCellDAObjArr[0].istep == 0)
			if (fabs(uebCellDAObjArr[0].currentModelDateTime - objEnKFArr.forecastDateTime) <= (0.5*ModelDt / 24))
			{
				cout << " Current datetime = " << CurrentDateTime[0] << " " << CurrentDateTime[1] << " " << CurrentDateTime[2] << " " << CurrentDateTime[3] << endl;
				cout << " Ens histYear start = " << objEnKFArr.ModelStartDateEns[0] << endl;
				cout << " Ens histYear end = " << objEnKFArr.ModelStartDateEns[0] + ensSize - 1 << endl;
			}

			for (int ensNum = 0; ensNum < ensSize; ensNum++)
			{
				boost::posix_time::ptime cDateTime(boost::gregorian::date(objEnKFArr.ModelStartDateEns[0] + ensNum, CurrentDateTime[1], CurrentDateTime[2]), //ptime(gregorian::date d, time_duration_type td)  //boost::gregorian::date(year_type y, month_type m, day_type d)
					boost::posix_time::hours(CurrentDateTime[3])); //  //boost::posix_time::ptime currentDateTime = currentDate +  get_property( g, graph_time_step_t() );	
				//cout << " CurrentDatetime = " << cDateTime << endl;
				//cout << " here -5 " << endl;
				vector< float > valT;
				for (int ifn = 0; ifn < 6; ifn++)
				{
					if (error = get_data_descript(g, ueb_forcNames[ifn], tempDesc)) return error;  //get_data_descript(GraphType& g, const char* name, ModelDataDescript& des)

					for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
					{//for indItr
						gd1dEns.readAtTime(cDateTime, inLoop[ifn], *indItr, valT);			//gd1d.readAtTime(curtime, des.name, *indItr, val);
						error = put_independent_basin_values(g, indItr->c_str(), valT, tempDesc);
						if (error)
						{//if
							throw std::runtime_error("put_independent_basin_values Error!");
						}//if
					}//for indItr 
					 //for (int ix = 0; ix < valT.size(); ix++) std::cout << " " << valT[ix];
				}
				if (error = get_values(g, ueb_forcNames, 6, &uebForcEnsArray[0])) return error;
				/*for (int ix = 0; ix < npix; ix++)
				std::cout << " " << uebForcEnsArray[npix + ix];
				std::cout << std::endl;*/

				if (useRDHMUnits)
				{
					//unit conversion  v (mm) during dT(min) :  m/hr => v * 60/(dt * 1000.0)	
					std::transform(uebForcEnsArray.begin(), uebForcEnsArray.begin() + npix, uebForcEnsArray.begin(), std::bind1st(std::multiplies<float>(), precConvFactor));
					//Tc = 5/9(Tf - 32) [](float inV) {return 5.0*(inV - 32.0) / 9.0;}
					std::transform(uebForcEnsArray.begin() + npix, uebForcEnsArray.begin() + 4 * npix, uebForcEnsArray.begin() + npix, [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
					//std::transform(uebTaminArray.begin(), uebTaminArray.end(), uebTaminArray.begin(), [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
					//std::transform(uebTamaxArray.begin(), uebTamaxArray.end(), uebTamaxArray.begin(), [](float inV) {return 5.0*(inV - 32.0) / 9.0; });
				}

				// 3.20.18 copy ensemble hist forc  copyEnsForcingHist(float *uebForcEnsArray, int ensNum, int cellRank, int numCells) 
				for (irank = 0; irank < npix; irank++)
					uebCellDAObjArr[irank].copyEnsForcingHist(&uebForcEnsArray[0], ensNum, irank, npix);
				/*=============== 9.3.18 for forecast period use grid closest to the obs point ================*/
				for (irank = 0; irank < nOpts; irank++)
				{
					daUebCellObjArr[irank].copyEnsForcingHist(&uebForcEnsArray[0], ensNum, objEnKFArr.Hc_hgVector[irank], npix);  //find corr grid
					//11.23.18 each obs points da without its observation (Leave one approach)
					daUebCellObjArr[irank + nOpts].copyEnsForcingHist(&uebForcEnsArray[0], ensNum, objEnKFArr.Hc_hgVector[irank], npix);  //find corr grid
				}
			}
			//put the actual current value--baseline year
			boost::posix_time::ptime cDateTime(boost::gregorian::date(CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2]), //ptime(gregorian::date d, time_duration_type td)  //boost::gregorian::date(year_type y, month_type m, day_type d)
				boost::posix_time::hours(CurrentDateTime[3])); //  //boost::posix_time::ptime currentDateTime = currentDate +  get_property( g, graph_time_step_t() );	
	//cout << " CurrentDatetime = " << cDateTime << endl;
			//cout << " here -5 " << endl;
			vector< float > valT;
			for (int ifn = 0; ifn < 6; ifn++)
			{
				if (error = get_data_descript(g, ueb_forcNames[ifn], tempDesc)) return error;  //get_data_descript(GraphType& g, const char* name, ModelDataDescript& des)

				for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
				{//for indItr
					gd1dDetr.readAtTime(cDateTime, inLoop[ifn], *indItr, valT);			//gd1d.readAtTime(curtime, des.name, *indItr, val);
					error = put_independent_basin_values(g, indItr->c_str(), valT, tempDesc);
					if (error)
					{//if
						throw std::runtime_error("put_independent_basin_values Error!");
					}//if
				}//for indItr 
			}
			/*if (error = get_values(g, ueb_forcNames, 6, uebForcEnsArray)) return error;
			for (int ix = 0; ix < npix; ix++)
			std::cout << " " << uebForcEnsArray[npix + ix];
			std::cout << std::endl;*/
			//cout << " here -3 " << endl;
//delete array
			//delete[] uebForcEnsArray;
			//uebForcEnsArray = NULL;
//cout << " here -2 " << endl;	
/*=============== 9.3.18 for forecast period use grid closest to the obs point ================*/
			/*Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> multivarNormalDistSamplesForc(3 * nOpts, ensSize);
			multivarNormalDistSamplesForc = objEnKFArr.std_norm_dist_Forc_Default.samples(ensSize);
			std::ofstream stdnormalSamplestxt("ensTempForcing_err.txt");
			stdnormalSamplestxt << multivarNormalDistSamplesForc << std::endl;

			std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > multivarNormalDistSamplesForcTVRH(3,
				Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(npix + nOpts, ensSize));
			multivarNormalDistSamplesForcTVRH[0] = objEnKFArr.norm_dist_1Mean_Default.samples(ensSize);
			multivarNormalDistSamplesForcTVRH[1] = objEnKFArr.norm_dist_1Mean_Default.samples(ensSize);
			multivarNormalDistSamplesForcTVRH[2] = objEnKFArr.norm_dist_1Mean_Default.samples(ensSize);
			std::ofstream stdnormalSamplestxtV("ensTempForcing_errT.txt");
			stdnormalSamplestxtV << multivarNormalDistSamplesForcTVRH[0] << std::endl;
			for (irank = 0; irank < nOpts; irank++)
				daUebCellObjArr[irank].copyEnsForcingMultiplier(irank + npix, ensSize, npix + nOpts, multivarNormalDistSamplesForc, multivarNormalDistSamplesForcTVRH);
			*/
			//call device run function	//using ens forcing from historical forcing
			callUEBrunHostEnsemble_HistForcing(threadsPerBlock, npix, ensSize, uebCellDAObjArr, nOpts, daUebCellObjArr);		//, &uebForcArrayDaPoint[0]);  //,uebEnKFDA objEnKFArr   
		}
		//computeRun_Time += (MPI::Wtime() - intermStart_Time);
		//intermStart_Time = MPI::Wtime();
//if (objEnKFArr.uebDAsacrutpix7_debugout1 == 1)
	//cout << " finished device compute tasks" << endl;
		//end of gpu call

		//EnKF
		//stateupdate arrays
		Eigen::RowVectorXf stateOutputArr(ensSize);
		Eigen::RowVectorXf stateOutputUpdateArr(ensSize + 1);
		std::vector<Eigen::RowVectorXf> Xh_stateObsSpace(nOpts, Eigen::RowVectorXf(ensSize));
		//2.20.18  xmrg outputs---last one for ens mean
		//std::vector< std::vector<float> > xmrgEnsembleArray(objEnKFArr.es_enseSize + 1, std::vector<float>(npix));
		//std::vector< std::vector<float> > xmrgDAEnsembleArray(ensSize + 1, std::vector<float>(npix));

		//ens and da outputs for time step --put obs points last 
		//float ***ensoutArray = create3DArrayblock_Contiguous(numEnsout - 1, npix + nOpts, ensSize + 1);
		float ***daoutArray = create3DArrayblock_Contiguous(numDAout - 1, npix + 2 * nOpts, ensSize + 1);
		/*std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > P_stateCovBackground(npix);
		std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > P_stateCovBackground_Points(npix);    //10.11.17 at obs points
		std::cout << "Saved cov. matrices" << std::endl;*/
//cout << " here 0 " << endl;

/*========hist forc-----
		int iranko;
		if (objEnKFArr.uebDAsacrutpix7_debugout3 == 1)
			std::cout << " State in obs. space: " << std::endl;
		for (irank = 0; irank < nOpts; irank++)
		{
			iranko = objEnKFArr.Hc_hgVector[irank];
			for (int ie = 0; ie < objEnKFArr.es_enseSize; ie++)
				Xh_stateObsSpace[irank](ie) = uebCellDAObjArr[iranko].stateVS[objEnKFArr.stateIndex][ie];  //stateOutputArr[iranko](ie);
			//
			if (objEnKFArr.uebDAsacrutpix7_debugout3 == 1)
			{
				std::cout << Xh_stateObsSpace[irank] << std::endl;
			}
		}
*//*======== 	using forc perturbation	=*/
        for (irank = 0; irank < nOpts; irank++)
		{
//cout << " here 00 " ;
			//copy states
			for (int ie = 0; ie < ensSize; ie++)
			{
				Xh_stateObsSpace[irank](ie) = daUebCellObjArr[irank].stateVS[objEnKFArr.stateIndex][ie];
//cout <<endl<< " here 01 ";
            }
			if (objEnKFArr.uebDAsacrutpix7_debugout3 == 1)
			{
				std::cout << Xh_stateObsSpace[irank] << std::endl;
			}
		}
/*=====================*/
//cout <<endl<< " here 1 " << endl;		
		Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> ensObservationErr(nOpts, ensSize);
		// generate randobs error only when there is assimilation
		if (fabs(objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA - 1] - uebCellDAObjArr[0].currentModelDateTime) <= (0.5*ModelDt / 24))    //call only when there is assimilation data---else no data (o) will be printed 8.5.16  ---<= 0.5 DT assimilate to the closest time step
		{
			ensObservationErr = objEnKFArr.norm_dist_0Mean_Default.samples(ensSize);
			std::ofstream stdnormalSamplestxtO("ensObsErr.txt");
			stdnormalSamplestxtO << ensObservationErr.transpose() << std::endl;
		}
		for (irank = 0; irank < npix; irank++)
		{
			//copy states
			for (int ie = 0; ie < ensSize; ie++)
				stateOutputArr(ie) = uebCellDAObjArr[irank].stateVS[objEnKFArr.stateIndex][ie];
//cout << " here 2 ";	
			if (objEnKFArr.uebDAsacrutpix7_debugout3 == 1)
			{
				if (irank == 0)
				{
					std::cout << std::endl << " cell 0 SWE ensemble outputs for time step: " << ueb_curModelTimeStep << std::endl;
					std::cout << stateOutputArr << std::endl;
				}
			}
			//ens output: 9.15.17.
			/*for (int icout = 0; icout < numEnsout - 1; icout++)
			{
				for (int ie = 0; ie < ensSize + 1; ie++)
					ensoutArray[icout][irank][ie] = uebCellDAObjArr[irank].stateVS[objEnKFArr.daEnsArr[icout].outvarIndex][ie];
			}*/
			if (fabs(objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA - 1] - uebCellDAObjArr[0].currentModelDateTime) <= (0.5*ModelDt / 24))    //call only when there is assimilation data---else no data (o) will be printed 8.5.16  ---<= 0.5 DT assimilate to the closest time step
			{
				objEnKFArr.runEnKF(objEnKFArr.daRegArray[uebCellDAObjArr[0].startIndexDA - 1], Xh_stateObsSpace, ensObservationErr, stateOutputArr, stateOutputUpdateArr);
				if (objEnKFArr.uebDAsacrutpix7_debugout3 == 1)
				{
					std::cout << std::endl << " sample updated state ensemble for time step: " << ueb_curModelTimeStep << " grid cell: " << irank << std::endl;
					std::cout << stateOutputUpdateArr << std::endl;       //is = danumStates*irank*size + is = danumStates*0*size + is					
				}
				//update background state for next step		
				uebCellDAObjArr[irank].updateBackgroundStates(objEnKFArr.es_enseSize, objEnKFArr.stateIndex, stateOutputUpdateArr);
			}	
			//da output: 9.15.17 apart from the state in EnKF the effect for other variables (e.g., melt) only appears next time step.
			for (int icout = 0; icout < numDAout - 1; icout++)
			{
				for (int ie = 0; ie < objEnKFArr.es_enseSize + 1; ie++)
					daoutArray[icout][irank][ie] = uebCellDAObjArr[irank].stateVS[objEnKFArr.daOutArr[icout].outvarIndex][ie];
			}
//cout << " here 3 ";	
			//2.20.18 copy xmrg da ens for da ens for Q simulation 			
			for (int ie = 0; ie < ensSize + 1; ie++)
			{
				//xmrgDAEnsembleArray[ie][irank] = uebCellDAObjArr[irank].stateVS[8][ie] * ModelDt;    //unit mm/dt  * 1000.0 * ModelDt;    //unit mm/dt uebCellDAObjArr[irank].OutVarValues[25] * 1000.0 * ModelDt;    //unit mm
				for (int i7 = 0; i7 < objEnKFArr.es_pfSize_Fact; i7++)
					pxvEns[irank * particleSize + objEnKFArr.es_pfSize_Fact * ie + i7] = uebCellDAObjArr[irank].stateVS[8][ie] * ModelDt;
			}
//cout << " here 3b ";	
        }
//cout << " here 4 " << endl;	
		for (irank = 0; irank < nOpts; irank++) //at obs points
		{
			//ens output: 9.15.17.
			/*for (int icout = 0; icout < numEnsout - 1; icout++)
			{
				for (int ie = 0; ie < ensSize + 1; ie++)
					ensoutArray[icout][irank + npix][ie] = daUebCellObjArr[irank].stateVS[objEnKFArr.daEnsArr[icout].outvarIndex][ie];
			}*/
			if (fabs(objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA - 1] - uebCellDAObjArr[0].currentModelDateTime) <= (0.5*ModelDt / 24))    //call only when there is assimilation data---else no data (o) will be printed 8.5.16  ---<= 0.5 DT assimilate to the closest time step
			{
				//1.15.18 ens already run no need to repeat
				objEnKFArr.runEnKF(objEnKFArr.daRegArray[uebCellDAObjArr[0].startIndexDA - 1], Xh_stateObsSpace, ensObservationErr, Xh_stateObsSpace[irank], stateOutputUpdateArr);
				//9.6.18 only assimilate at grid cell obs
				//objEnKFArr.runEnKFopt(irank, objEnKFArr.daRegArray[uebCellDAObjArr[0].startIndexDA - 1], Xh_stateObsSpace[irank], ensObservationErr, Xh_stateObsSpace[irank], stateOutputUpdateArr);
				if (objEnKFArr.uebDAsacrutpix7_debugout3 == 1)
				{
					std::cout << " sample updated state ensemble for time step: " << ueb_curModelTimeStep << " observation point: " << irank << std::endl;
					std::cout << stateOutputUpdateArr << std::endl;       //is = danumStates*irank*size + is = danumStates*0*size + is					
				}
				//update background state for next step		
				daUebCellObjArr[irank].updateBackgroundStates(objEnKFArr.es_enseSize, objEnKFArr.stateIndex, stateOutputUpdateArr);
	
			}
//cout << " here 5 ";	
//da output 9.15.17 apart from the state in EnKF the effect for other variables (e.g., melt) only appears next time step.
			for (int icout = 0; icout < numDAout - 1; icout++)
			{
				for (int ie = 0; ie < objEnKFArr.es_enseSize + 1; ie++)
					daoutArray[icout][irank + npix][ie] = daUebCellObjArr[irank].stateVS[objEnKFArr.daOutArr[icout].outvarIndex][ie];
			}

			char pntoStr[256];
			sprintf(pntoStr, "%d", irank);
			char outObsPoint[256];
			strcpy(outObsPoint, "obsSNOTELpoint");    // one nc file stores arrays for all variables
			strcat(outObsPoint, pntoStr);
			strcat(outObsPoint, ".txt");
			//iranko = objEnKFArr.Hc_hgVector[irank];
			//uebCellDAObjArr[iranko].printPointOutputs((const char*)outObsPoint);
			daUebCellObjArr[irank].printPointOutputs((const char*)outObsPoint);
		}	
		//11.23.18 each obs points da without its observation (Leave one approach)
		for (irank = 0; irank < nOpts; irank++) //at obs points
		{
			if (fabs(objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA - 1] - uebCellDAObjArr[0].currentModelDateTime) <= (0.5*ModelDt / 24))    //call only when there is assimilation data---else no data (o) will be printed 8.5.16  ---<= 0.5 DT assimilate to the closest time step
			{
				objEnKFArr.runEnKF_LeaveOneOut(irank, objEnKFArr.daRegArray[uebCellDAObjArr[0].startIndexDA - 1], Xh_stateObsSpace, ensObservationErr, Xh_stateObsSpace[irank], stateOutputUpdateArr);
				//update background state for next step		
				daUebCellObjArr[irank + nOpts].updateBackgroundStates(objEnKFArr.es_enseSize, objEnKFArr.stateIndex, stateOutputUpdateArr);
			}
			//cout << " here 5 ";	
			for (int icout = 0; icout < numDAout - 1; icout++)
			{
				for (int ie = 0; ie < objEnKFArr.es_enseSize + 1; ie++)
					daoutArray[icout][nOpts + irank + npix][ie] = daUebCellObjArr[irank + nOpts].stateVS[objEnKFArr.daOutArr[icout].outvarIndex][ie];
			}
			char pntoStr[256];
			sprintf(pntoStr, "%d", irank);
			char outObsPoint[256];
			strcpy(outObsPoint, "obsSNOTELpoint");    // one nc file stores arrays for all variables
			strcat(outObsPoint, pntoStr);
			strcat(outObsPoint, "_leaveOut.txt");
			//iranko = objEnKFArr.Hc_hgVector[irank];
			//uebCellDAObjArr[iranko].printPointOutputs((const char*)outObsPoint);
			daUebCellObjArr[irank + nOpts].printPointOutputs((const char*)outObsPoint);
		}
//cout << " here 6 " << endl;		
		if (fabs(objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA - 1] - uebCellDAObjArr[0].currentModelDateTime) <= (0.5*ModelDt / 24))    //call only when there is assimilation data---else no data (o) will be printed 8.5.16  ---<= 0.5 DT assimilate to the closest time step
		{
			//std::cout << " done DA update this da time: " << objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA - 1] <<
			//	" curr model time: " << uebCellDAObjArr[0].currentModelDateTime << std::endl;
			//std::cout << " next DA time: " << objEnKFArr.daTcorrArr[uebCellDAObjArr[0].startIndexDA] << std::endl;
			uebCellDAObjArr[0].updateDaArray = true;            //copy the next observation 			
		}   
//#3.12.18 for tracking sac and rx states during ensemble run	
		/*int ensRemainder = 0;
		std::vector<float> ensXmrgMultiplier(7, 1.0);  // = { 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.0 };
		std::vector<float> currXmrgArray(npix);*/

		//float *QoutEnsArr = new float[particleSize];
		float *QoutDAArr = new float[particleSize];
		std::vector<int> rx7ndxArray(npix);
		std::vector< std::vector<float> > rx7DepthArray(npix, std::vector<float>(particleSize));		// objEnKFArr.es_enseSize + 1));
		//std::vector< Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > sacrealStateArray(npix, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(6, particleSize));		// objEnKFArr.es_enseSize + 1));
		std::vector < Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > rx7AreaCiArray(npix, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(4, particleSize));		// objEnKFArr.es_enseSize + 1));
		if (error = get_value(g, "ndx", &rx7ndxArray[0])) return error;
		if (error = get_value(g, "rx7DepthArray", &rx7DepthArray[0])) return error;
		//if (error = get_value(g, "sacrealStateArray", &sacrealStateArray[0])) return error;
		if (error = get_value(g, "rx7AreaCiArray", &rx7AreaCiArray[0])) return error;
		
//1.9.19 EDmnd
		int year, month, day, hour;
		int dtm;
		if (error = get_year_month_day_hour(g, year, month, day, hour)) return error;
		if (error = get_time_step_minutes(g, dtm)) return error;
		bt::shared_array< float > pe;
		bt::shared_array< float > pe_adj;
		retrieve_value(g, "pPE", pe);
		retrieve_value(g, "pPEAdj", pe_adj);
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
			sacGridArray[pix].aesc = 0;

			pix++;
		}
		
		do_sac_EnsHost_Ens(threadsPerBlock, npix, particleSize, sacGridArray, pxvEns);
//cout << " here 7 " << endl;

		for (int ie = 0; ie < particleSize; ie++)		// objEnKFArr.es_enseSize + 1; ie++)
		{
			//####update prev sac and rx7 states			
			irank = 0;
			for (itr = selectedPixels.begin(); itr != selectedPixels.end(); ++itr)
			{//for itr for (irank = 0; irank < npix; irank++)	{
				for (int indx = 0; indx < rx7ndxArray[irank]; indx++)
				{
					//std::wcout << " here 9" << std::endl;
					namePos = nameIndex.find("areac" + boost::lexical_cast<string>(indx + 1));
					if (namePos == nameIndex.end())
					{
						std::cout << " Cannopt put areac val: Property areac" << indx + 1 << " not found at grid cell " << irank << std::endl;
						//std::getchar();
						return PROPERTY_NOT_FOUND;
					}//if
					else
					{
						//std::wcout << " here 11" << std::endl;
						attMap[*itr][namePos->second] = rx7AreaCiArray[irank](indx, ie);
					}
					//dx_st_g[indx*npix + irank] = rx7AreaCiArray[irank](indx, ie);
				}//for indx
				//std::wcout << " here 9" << std::endl;
				//if (error = put_value(g, &rx7DepthArray[ie][0], "depth")) return error;
				namePos = nameIndex.find("depth");
				if (namePos == nameIndex.end())
				{
					std::cout << " Cannot put var depth at grid cell " << irank << std::endl;
					//std::getchar();
					return PROPERTY_NOT_FOUND;
				}//if
				else
				{
					//std::wcout << " here 11" << std::endl;
					attMap[*itr][namePos->second] = rx7DepthArray[irank][ie];
				}
				//rut_st_g[irank] = rx7DepthArray[irank][ie];
				/*for (int irstn = 0; irstn < 6; irstn++)
				{
					//std::wcout << " here 9" << std::endl;
					namePos = nameIndex.find(sac_real_stnames[irstn]);
					if (namePos == nameIndex.end())
					{
						std::cout << " Cannot put " << sac_real_stnames[irstn] << " at grid cell " << irank << std::endl;
						//std::getchar();
						return PROPERTY_NOT_FOUND;
					}//if
					else
					{
						//std::wcout << " here 11" << std::endl;
						attMap[*itr][namePos->second] = sacrealStateArray[irank](irstn, ie);
						//sac_st_real_Ens[irstn * npix + irank] = sacrealStateArray[irank](irstn, ie);
					}
				}//for 	*/
				//###-----update sac output
				attMap[*itr][nameIndex["surfaceFlow"]] = boost::any(sacGridArray[irank].surfEns[ie]);
				attMap[*itr][nameIndex["subsurfaceFlow"]] = boost::any(sacGridArray[irank].grndEns[ie]);
				/*namePos = nameIndex.find("surfaceFlow");
				if (namePos == nameIndex.end())
				{
					std::cout << " Cannot put var surfaceFlow at grid cell " << irank << std::endl;
					//std::getchar();
					return PROPERTY_NOT_FOUND;
				}//if
				else
				{
					//std::wcout << " here 11" << std::endl;
					attMap[*itr][namePos->second] = sacGridArray[irank].surfEns[ie];
				}
				//std::wcout << " here 10" << std::endl;
				namePos = nameIndex.find("subsurfaceFlow");
				if (namePos == nameIndex.end())
				{
					std::cout << " Cannot put var subsurfaceFlow at grid cell " << irank << std::endl;
					//std::getchar();
					return PROPERTY_NOT_FOUND;
				}//if
				else
				{
					//std::wcout << " here 11" << std::endl;
					attMap[*itr][namePos->second] = sacGridArray[irank].grndEns[ie];
				}*/

				irank++;
			}
			//update xmrg
			//if (error = put_value(g, &xmrgDAEnsembleArray[ie / 7][0], "xmrg")) return error;
			//if (error = put_value(g, &xmrgDAEnsembleArray[ie / 7][0], "rmlt")) return error;

/*cout << endl << " Before doSac: Subsurface flow " << endl;*/
			//if (error = sacFuncInside(g)) return error;
			//if (error = sacFuncInside_Ens(g, threadsPerBlock, xmrgDAEnsembleArray[ie / 7], sacGridArray, 0)) return error;


/*cout << endl << " After doSac: Subsurface flow " << endl;
*/
			if (error = rutpixFuncInside(g)) return error;
			/*
						if (error = calrutpixFuncInside_Ens(g, surfaceFlowEns, subsurfaceFlowEns)) return error;
						do_route_(&surflow_g[0], &subflow_g[0], &qx_g[0],
							&dtsec_g,
							&rut_par_g[0], &npix, &pixarea_g[0], &ndx_g[0], &length_g[0],
							&idown_g[0], &rut_st_g[0], &dx_st_g[0], &num_fpt_g[0], &nordi_g[0],
							&upstream_replace_g[0], &upqx_g[0], &chan_loss_g[0]
							);
						if (error = calrutpixFuncAfter_Ens(g)) return error;
			*/
			//std::cout << " Done Rutpix7 run " << endl;

						//### 2.20.18 till better way found
			std::vector<float> qVect;
			if (error = getQ(g, qVect)) return error;
			//int nfpt; if (error = get_nfpt(g, nfpt)) return error;
			//for (int ib = 0; ib < nfpt; ib++) QoutDAArr[ie] =  qVect[ib]; 
			//QoutEnsArr[ie] = qVect[0];
			/*if (qVect[0] < 0 || isnan(qVect[0]) || isinf(qVect[0]))
			{
				QoutDAArr[ie] = 0.0;
				qVect[0] = 0.0;
				if (error = putQ(g, qVect)) return error;
			}
			else
			{*/
			QoutDAArr[ie] = qVect[0];
			//std::cout << " Q= "<< qVect[0] << " ";
			//if (uebCellDAObjArr[0].istep == 0) std::cout << "number of outlet points: " << qVect.size() << std::endl; 			
			//#### Save ens sac and rx7 states for next step
			irank = 0;
			for (itr = selectedPixels.begin(); itr != selectedPixels.end(); ++itr)
			{//for itr for (irank = 0; irank < npix; irank++)	{
				for (int indx = 0; indx < rx7ndxArray[irank]; indx++)
				{
					//std::wcout << " here 9" << std::endl;
					namePos = nameIndex.find("areac" + boost::lexical_cast<string>(indx + 1));
					if (namePos == nameIndex.end())
					{
						std::cout << " Cannopt save areac val: Property areac" << indx + 1 << " not found at grid cell " << irank << std::endl;
						//std::getchar();
						return PROPERTY_NOT_FOUND;
					}//if
					else
					{
						rx7AreaCiArray[irank](indx, ie) = boost::any_cast<float>(attMap[*itr][namePos->second]);
					}
					//rx7AreaCiArray[irank](indx, ie) = dx_st_g[indx*npix + irank];
				}//for indx
//std::cout << " here 9_1 " ;
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
					rx7DepthArray[irank][ie] = boost::any_cast<float>(attMap[*itr][namePos->second]);
				}
				irank++;
			}
			//cout << " here 10 " << endl;
		}
//cout << endl << " here 13 " << endl;	

		if (fabs(objEnKFArr.daTcorrArrQ[uebCellDAObjArr[0].startIndexDAQ - 1] - uebCellDAObjArr[0].currentModelDateTime) <= (0.5*ModelDt / 24))    //call only when there is assimilation data---else no data (o) will be printed 8.5.16  ---<= 0.5 DT assimilate to the closest time step
		{
			
//##4.9.18:	//particle weight based on norm dist
			//Eigen::RowVectorXf partWeightRx7(particleSize);
			std::vector< Eigen::RowVectorXf> partWeight(npix, Eigen::RowVectorXf(particleSize));
			if (error = get_value(g,  "usrParticleWeights", &partWeight[0])) return error;

			for (int ie = 0; ie < particleSize; ie++)
			{
				float qDiff = QoutDAArr[ie] - objEnKFArr.daQstrArray(uebCellDAObjArr[0].startIndexDAQ - 1);
				partWeight[0](ie) *= exp(-0.5 * qDiff * qDiff / (objEnKFArr.obsQErrStdev * objEnKFArr.obsQErrStdev));
				//for rutpix7
				//partWeightRx7(ie) = exp(-0.5 * qDiff * qDiff / (objEnKFArr.obsQErrStdev * objEnKFArr.obsQErrStdev));
			}
			//cout << " Weights from lik func " << std::endl << partWeight << std::endl;
			float sumWeights = partWeight[0].sum();
			//normalize and cummWeight
			partWeight[0] = partWeight[0] / sumWeights;
			//cout << " Normalied Weights " << std::endl << partWeight[0] << std::endl;
//#7.30.18 resample/update
			if ((CurrentDateTime[2] % objEnKFArr.qUpdateFreq) == 0 )		//CurrentDateTime[2] == 1 || CurrentDateTime[2] == 7 || CurrentDateTime[2] == 15 || CurrentDateTime[2] == 22)   //(fabs(uebCellDAObjArr[0].currentModelDateTime - objEnKFArr.forecastDateTime) <= (0.5*ModelDt / 24))
			{
				//for sac update
				std::vector< Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > sacrealStateArray_update(npix, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(6, particleSize));		// objEnKFArr.es_enseSize + 1));
				//for q update
				std::vector< std::vector<float> > rx7DepthArray_update(npix, std::vector<float>(particleSize));		// objEnKFArr.es_enseSize + 1));			
				std::vector < Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > rx7AreaCiArray_update(npix, Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(4, particleSize));		// objEnKFArr.es_enseSize + 1));
				float *QoutDAArr_update = new float[particleSize];
				
				Eigen::RowVectorXf cummWeight(particleSize);
				cummWeight(0) = partWeight[0](0);
				float sumSqWeights = 0.0;
				for (int ie = 1; ie < particleSize; ie++)
				{
					sumSqWeights += (partWeight[0](ie) * partWeight[0](ie));
					//partWeight[ie] = partWeight[ie] / sumWeights;
					cummWeight(ie) = cummWeight(ie - 1) + partWeight[0](ie);
				}
				int effectiveParticleSize = int(1.0 / sumSqWeights);
				cout << " Effective particle size " << effectiveParticleSize << std::endl;
				//if(cummWeight(particleSize - 1) > 1.0)
				cummWeight(particleSize - 1) = 1.0;
				cout << " Cummulative Weights " << cummWeight << std::endl;
				/*******************************************************************************************************************************/
				/*****ToDO: Needs revision: ####systematic resampling:
				*/
				std::vector<int> weightIndices(particleSize, particleSize - 1);
				int ipx = 0, jpx = 0;
				while ((ipx < particleSize) && (jpx < particleSize))
				{
					if (((float)ipx / particleSize) < cummWeight(jpx))
					{
						weightIndices[ipx] = jpx;
						ipx++;
					}
					else
					{
						jpx++;
					}
					if (jpx >= particleSize)
						cout << " Warning index = " << jpx << " exceeded max = " << particleSize << std::endl;
				}
				/******************************************************************************************************************************/
				cout << " Selected Weight indices " << weightIndices << std::endl;
				for (int ie = 0; ie < particleSize; ie++)
				{
					QoutDAArr_update[ie] = QoutDAArr[weightIndices[ie]];
					irank = 0;
					for (itr = selectedPixels.begin(); itr != selectedPixels.end(); ++itr)
					{//for itr for (irank = 0; irank < npix; irank++)	{
						for (int indx = 0; indx < rx7ndxArray[irank]; indx++)
						{
							rx7AreaCiArray_update[irank](indx, ie) = rx7AreaCiArray[irank](indx, weightIndices[ie]);
						}//for indx
						
						for (int irstn = 0; irstn < 6; irstn++)
						{
							//sacrealStateArray_update[irank](irstn, ie) = sacrealStateArray[irank](irstn, weightIndices[ie]);
							sacrealStateArray_update[irank](irstn, ie) = sacGridArray[irank].sacstEns[weightIndices[ie]][irstn];
						}//for
						rx7DepthArray_update[irank][ie] = rx7DepthArray[irank][weightIndices[ie]];
						irank++;
					}
				}
				// update SAC-SMA
                //8.4.18 ==========use small perturb to avoid sample impoverishment =================================================
				std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > multivarNormalDistSamplesSACSMA(6,
					Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>(npix, particleSize));
				for (int iforc = 0; iforc < 6; iforc++)
				{
					multivarNormalDistSamplesSACSMA[iforc] = objEnKFArr.norm_dist_SACRX.samples(particleSize);
				}
				std::ofstream stdnormalSamplestxtSACSMA("particles_perturb_mult.txt");
				stdnormalSamplestxtSACSMA << multivarNormalDistSamplesSACSMA[0] << std::endl;
				//
				for (irank = 0; irank < npix; irank++)
				{
					for (int ie = 0; ie < particleSize; ie++)
					{
						for (int irstn = 0; irstn < 6; irstn++)
						{
							//sacrealStateArray_update[irank](irstn, ie) = sacrealStateArray[irank](irstn, weightIndices[ie]);
							sacGridArray[irank].sacstEns[ie][irstn] = sacrealStateArray_update[irank](irstn, ie) * fabs(multivarNormalDistSamplesSACSMA[irstn](irank, ie)); //8118 mult.by pert.fact;
						}//for
					}
				}
				//cout << " here 14.b " << std::endl;

				//update UEB

				//cout << " here 14.c " << std::endl;
				if (error = put_value(g, &rx7DepthArray_update[0], "rx7DepthArray")) return error;
				//if (error = put_value(g, &sacrealStateArray_update[0], "sacrealStateArray")) return error;
				if (error = put_value(g, &rx7AreaCiArray_update[0], "rx7AreaCiArray")) return error;
				//q
				retvalue = Write1DVector_to2DNC((const char*)objEnKFArr.daOutArr[numDAout - 1].outfName, (const char*)objEnKFArr.daOutArr[numDAout - 1].symbol, ueb_curModelTimeStep, particleSize, QoutDAArr_update);          //, worldComm, worldInfo);	
																																																								//																																																				//
				delete[] QoutDAArr_update;
				QoutDAArr_update = NULL;
				for (int ie = 0; ie < particleSize; ie++)		// objEnKF.es_enseSize + 1; ie++)
				{
					//initialize particle weights
					partWeight[0](ie) = 1.0;
				}
				//if (error = put_value(g, &partWeight[0], "usrParticleWeights")) return error;
				/*std::cout << std::setprecision(15) << " done SAC-SMA <- Q update this time: " << objEnKFArr.daTcorrArrQ[uebCellDAObjArr[0].startIndexDAQ - 1] <<
				" model time: " << uebCellDAObjArr[0].currentModelDateTime 
				 << " Current datetime = " << CurrentDateTime[0] << " " << CurrentDateTime[1] << " " << CurrentDateTime[2] << " " << CurrentDateTime[3] << endl;
				*/
				//std::cout << std::setprecision(15) << " next Qstr DA time step: " << uebCellDAObjArr[0].startIndexDAQ << " next daQ time: " << objEnKFArr.daTcorrArrQ[uebCellDAObjArr[0].startIndexDAQ] << std::endl;				
			}
			else
			{
				//if (error = put_value(g, &partWeight[0], "usrParticleWeights")) return error;
				if (error = put_value(g, &rx7DepthArray[0], "rx7DepthArray")) return error;
				//if (error = put_value(g, &sacrealStateArray[0], "sacrealStateArray")) return error;
				if (error = put_value(g, &rx7AreaCiArray[0], "rx7AreaCiArray")) return error;
				//save Q
				retvalue = Write1DVector_to2DNC((const char*)objEnKFArr.daOutArr[numDAout - 1].outfName, (const char*)objEnKFArr.daOutArr[numDAout - 1].symbol, ueb_curModelTimeStep, particleSize, QoutDAArr);          //, worldComm, worldInfo);
			}
			//
			if (error = put_value(g, &partWeight[0], "usrParticleWeights")) return error;			
		    //std::cout << std::setprecision(15) << " done Qstr DA update this time: " << objEnKFArr.daTcorrArrQ[uebCellDAObjArr[0].startIndexDAQ - 1] <<
            //	" model time: " << uebCellDAObjArr[0].currentModelDateTime << std::endl;
			//std::cout << std::setprecision(15) << " next Qstr DA time step: " << uebCellDAObjArr[0].startIndexDAQ << " next daQ time: " << objEnKFArr.daTcorrArrQ[uebCellDAObjArr[0].startIndexDAQ] << std::endl;
			uebCellDAObjArr[0].updateDaArrayQ = true;            //copy the next observation         
        }
		else
		{
			if (error = put_value(g, &rx7DepthArray[0], "rx7DepthArray")) return error;
			//if (error = put_value(g, &sacrealStateArray[0], "sacrealStateArray")) return error;
			if (error = put_value(g, &rx7AreaCiArray[0], "rx7AreaCiArray")) return error;
			//save Q
			retvalue = Write1DVector_to2DNC((const char*)objEnKFArr.daOutArr[numDAout - 1].outfName, (const char*)objEnKFArr.daOutArr[numDAout - 1].symbol, ueb_curModelTimeStep, particleSize, QoutDAArr);          //, worldComm, worldInfo);
		}
//cout << " here 15a " << endl;
		//put ens da outputs to nc
//cout << " here 15.b " << endl;
		/*for (int icout = 0; icout < numEnsout - 1; icout++)
		{
			retvalue = Write2DSlub_to3DNC((const char*)objEnKFArr.daEnsArr[icout].outfName, (const char*)objEnKFArr.daEnsArr[icout].symbol, ueb_curModelTimeStep, nOpts + npix, ensSize + 1, ensoutArray[icout]);   //, worldComm, worldInfo);
		}*/
		for (int icout = 0; icout < numDAout - 1; icout++)
		{
			retvalue = Write2DSlub_to3DNC((const char*)objEnKFArr.daOutArr[icout].outfName, (const char*)objEnKFArr.daOutArr[icout].symbol, ueb_curModelTimeStep, 2 * nOpts + npix, ensSize + 1, daoutArray[icout]);   //, worldComm, worldInfo);
		}
		//retvalue = Write1DVector_to2DNC((const char*)objEnKFArr.daEnsArr[numEnsout - 1].outfName, (const char*)objEnKFArr.daEnsArr[numEnsout - 1].symbol, ueb_curModelTimeStep, particleSize, QoutEnsArr);          //, worldComm, worldInfo);
//cout << " here 15 " << endl;
		for (irank = 0; irank < 2 * objEnKFArr.numObsPoints; irank++)
			daUebCellObjArr[irank].updateSimTime();
		if (error = put_value(g, daUebCellObjArr, "daUebCellObjArray"))	return error;	
//std::cout << " saved updated da objects and for time step " << uebCellDAObjArr[0].istep << std::endl;
		delete3DArrayblock_Contiguous(daoutArray);
		daoutArray = NULL;
		//delete3DArrayblock_Contiguous(ensoutArray);
		//ensoutArray = NULL;
//cout << " here 15b " << endl;
		//delete[] QoutEnsArr;
		//QoutEnsArr = NULL;
		delete[] QoutDAArr;
		QoutDAArr = NULL;
//cout << " here 15.c " << endl;
	} //if (uebCellDAObjArr[0].daAssimlate)
	else
	{	//run ueb--for cases where there is no da
		/*for (irank = 0; irank < npix; irank++)
		{
			uebCellDAObjArr[irank].setForcingAndRadiationParamterization(uebForcArray[irank], uebForcArray[npix + irank], uebForcArray[2*npix + irank],
			uebForcArray[3*npix + irank], uebForcArray[4*npix + irank], uebForcArray[5*npix + irank]);
			//uebCellDAObjArr[irank].runUEB();
		}*/
//cout <<endl<< " here 16a ";
		std::vector<float>  xmrgDAEnsembleArray(npix);
		//gpu control	
		//threadsPerBlock = 128;
		callUEBrunHost(threadsPerBlock, npix, uebCellDAObjArr, &uebForcArray[0]);

//window outputs: note expect small window size, so write all
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
				uebCellDAObjArr[irank].printPointOutputs((const char*)outObsPoint);
			}
		}		
//cout <<endl<< " here 16a2 ";		
		for (irank = 0; irank < npix; irank++)
		{
			xmrgDAEnsembleArray[irank] = uebCellDAObjArr[irank].OutVarValues[25] * ModelDt;    //unit mm/dt  * 1000.0 * ModelDt;    //unit mm/dt uebCellDAObjArr[irank].OutVarValues[25] * 1000.0 * ModelDt;    //unit mm
		}
//cout << " here 16b ";	
//*----------------------------- TODO: write point outputs
/*		for (irank = 0; irank < nOpts; irank++) //at obs points
		{
			char pntoStr[256];
			sprintf(pntoStr, "%d", irank);
			char outObsPoint[256];
			strcpy(outObsPoint, "obsSNOTELpoint");    // one nc file stores arrays for all variables
			strcat(outObsPoint, pntoStr);
			strcat(outObsPoint, ".txt");
			daUebCellObjArr[irank].printPointOutputs((const char*)outObsPoint);
		}
------------------------------------------*/
		//update xmrg
		if (error = put_value(g, &xmrgDAEnsembleArray[0], "xmrg")) return error;
		if (error = put_value(g, &xmrgDAEnsembleArray[0], "rmlt")) return error;
//cout << " here 16c ";
       
        //threadsPerBlock = 255;
		//if (error = sacFuncInside(g)) return error;	
		if (error = sacFuncInside_Ens(g, threadsPerBlock, xmrgDAEnsembleArray, sacGridArray, 0)) return error;
//std::cout << "  Done SAC run "; // << std::endl;

		if (error = rutpixFuncInside(g)) return error;
//std::cout << " Done Rutpix7 run "; // << endl;

		/*if (error = calsacFuncInside(g)) return error;
		if (error = calsacFuncAfter(g)) return error;
std::cout << "  Done SAC run "; // << std::endl;
		if (error = calrutpixFuncInside(g)) return error;
		if (error = calrutpixFuncAfter(g)) return error;		
std::cout << " Done Rutpix7 run "; // << endl;*/

		//### 2.20.18 till better way found
		float *QoutDAArr1 = new float[1];
		std::vector<float> qVect;
		if (error = getQ(g, qVect)) return error;

		QoutDAArr1[0] = qVect[0];
		retvalue = Write1DVector_to2DNC((const char*)objEnKFArr.daOutArr[numDAout - 1].outfName, (const char*)objEnKFArr.daOutArr[numDAout - 1].symbol, ueb_curModelTimeStep, 1, QoutDAArr1);          //, worldComm, worldInfo);
																																																	   //
		delete[] QoutDAArr1;
		QoutDAArr1 = NULL;

		//int nfpt; if (error = get_nfpt(g, nfpt)) return error;
		//for (int ib = 0; ib < nfpt; ib++) QoutDAArr[ie] =  qVect[ib]; 
//cout << " here 16.d " << endl;
	}
//cout <<" here 17 " << endl;		
	//save U, W, rmlt 
	for (irank = 0; irank < npix; irank++)
	{
		UbArray[irank] = uebCellDAObjArr[irank].OutVarValues[16];           // KJ/m^2
		WArray[irank] = uebCellDAObjArr[irank].OutVarValues[17];				//unit mm    //* 1000.0
		rmltArray[irank] = uebCellDAObjArr[irank].OutVarValues[25] * ModelDt;    //unit mm/ dt   * 1000.0 * ModelDt;    //unit mm
		xmrgArray[irank] = uebCellDAObjArr[irank].OutVarValues[25] * ModelDt;    //unit mm/ dt  *1000.0 * ModelDt;    //unit mm
	}
//std::cout << "here 17b"<<std::endl;
	//---UEB outputs: 11.30.17 This keeps compatibility to UEB C++/Fortran versions
	for (iuout = 0; iuout < 71; iuout++)
	{
		//uebOutputArray[iuout][irank] = uebCellDAObjArr[irank].OutVarValues[iuout];
		for (irank = 0; irank < npix; irank++)
			uebOutputArray[iuout * npix + irank] = uebCellDAObjArr[irank].OutVarValues[iuout];
//std::cout << "here 18 ";
	}
//std::cout << std::endl << " here 19 "<< std::endl;
	//UEB outputs: This keeps compatibility to UEB C++/Fortran versions
	const char* uebVars[] = { "ueb_Year", "ueb_Month", "ueb_Day", "ueb_dHour", "ueb_atff", "ueb_HRI", "ueb_Eacl", "ueb_Ema", "ueb_cosZen", "ueb_Ta", "ueb_P", "ueb_V", "ueb_RH",
		"ueb_Qsi", "ueb_Qli", "ueb_Qnet", "ueb_Us", "ueb_SWE", "ueb_tausn", "ueb_Pr", "ueb_Ps", "ueb_Alb", "ueb_QHs", "ueb_QEs", "ueb_Es", "ueb_SWIT", "ueb_QMs", "ueb_Q",
		"ueb_FM", "ueb_Tave", "ueb_TSURFs", "ueb_cump", "ueb_cumes", "ueb_cumMr", "ueb_NetRads", "ueb_smelt", "ueb_refDepth", "ueb_totalRefDepth", "ueb_cf", "ueb_Taufb",
		"ueb_Taufd", "ueb_Qsib", "ueb_Qsid", "ueb_Taub", "ueb_Taud", "ueb_Qsns", "ueb_Qsnc", "ueb_Qlns", "ueb_Qlnc", "ueb_Vz", "ueb_Rkinsc", "ueb_Rkinc", "ueb_Inmax",
		"ueb_intc", "ueb_ieff", "ueb_Ur", "ueb_Wc", "ueb_Tc", "ueb_Tac", "ueb_QHc", "ueb_QEc", "ueb_Ec", "ueb_Qpc", "ueb_Qmc", "ueb_Mc", "ueb_FMc", "ueb_SWIGM",
		"ueb_SWISM", "ueb_SWIR", "ueb_errMB", "ueb_Trange" };
	/*for (iuout = 0; iuout < 71; iuout++)
	{
		if (error = put_value(g, &uebOutputArray[iuout][0], uebVars[iuout])) return error;
	}*/
	if (error = put_values(g, 71, &uebOutputArray[0], uebVars)) return error;
//std::cout << "here 20 " << std::endl;
	//put outputs (that may be used by other models, e.g., xmrg by sacsma)
	if (error = put_value(g, &xmrgArray[0], "xmrg")) return error;
	if (error = put_value(g, &rmltArray[0], "rmlt")) return error;
	if (error = put_value(g, &WArray[0], "uebW")) return error;
	if (error = put_value(g, &UbArray[0], "uebUb")) return error;
	//attMap[*itr][nameIndex["rmlt"]] = boost::any(rmlt);
//std::cout << "here 21 " << std::endl;

	//set next time step
	for (irank = 0; irank < npix; irank++)	
		uebCellDAObjArr[irank].updateSimTime();	
	//save ueb obj to g (pixgraph) obj
	if (error = put_value(g, &uebCellDAObjArr[0], "uebCellDAObjArray")) return error;
	//if (error = put_value(g, &objEnKFArr, "uebEnKFDAObjArr")) return error;
	if (error = save_value(g, "uebEnKFDAObj", objEnKFArr)) return error;
	//
	if (error = put_value(g, &sacGridArray[0], "des_sacGridArray")) return error;
	//std::cout << " saved updated ueb objects and outputs for time step " << uebCellDAObjArr[0].istep << std::endl;
//std::cout << "here 22 " << std::endl;
	
	int prog_istep = uebCellDAObjArr[0].istep;
	int prog_numTotalTimeSteps = uebCellDAObjArr[0].numTotalTimeSteps;
	//delete point forc array
	//delete[] uebForcArray;
	//uebForcArray = NULL;
	delete[] uebCellDAObjArr;
	uebCellDAObjArr = NULL;
//std::cout << "here 23 " << std::endl;
	//delete[] uebForcArrayDaPoint;
	//uebForcArrayDaPoint = NULL;	
	delete[] daUebCellObjArr;
	daUebCellObjArr = NULL;

	delete[] sacGridArray;
	sacGridArray = NULL;
//std::cout << "here 24 " << std::endl;
	delete[] pxvEns;
	pxvEns = NULL;

	//std::cout << " freed ens/da ooutput arrays" << std::endl;

	//progress is calculated and written here		
	//if (rank == 0) 
	if ((prog_istep % 120) == 0)
		std::cout << " Percent completed: " << ((float)prog_istep / prog_numTotalTimeSteps)*100.0 << " %" << std::endl;
	
	//	cerr << "leaving uebFuncInsideLoop" << endl;
	return OK;

}//uebDAsacrutpix7FuncInsideLoop
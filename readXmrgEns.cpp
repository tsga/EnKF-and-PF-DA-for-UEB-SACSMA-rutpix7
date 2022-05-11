#include <vector>
#include <iostream>
#include <sstream>
#include "boost/lexical_cast.hpp"
#include "boost/shared_array.hpp"
#include "getSelectedVerticesInConnOrder.hpp"
#include "hlrms_interfaces.hpp"
#include "getFromDeck.hpp"
#include "makeModelDataDescription.h"

#include <functional>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/format.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/fstream.hpp"
#include "PixelGraph.hpp"
#include "ModelBase.hpp"
#include "hlrms_interfaces.hpp"
#include "getInputDeck.h"
#include "accumulateGrids.h"
#include "selectProperties.h"
#include "GridData1D.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "GridData1D.h"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"


	using namespace Ohd_Hydro;
	using namespace std;

	extern const Ohd_Hydro::DECK deck;

	namespace bt = boost;
	namespace gg = bt::gregorian;
	namespace pt = bt::posix_time;
	namespace fs = bt::filesystem;
	//
	PixelGraph& allPixels;
	ModelDataDescript const& des;
	//
	// Get 1D XMRGs
	//
	/*vector< string > inLoop;
	ModelPropertyMap::value_type::second_type inLoopDes;
	ModelPropertyMap::value_type::second_type::iterator desItr;
	selectProperties(deck, mptInsideLoop, inLoopDes);
	for (desItr = inLoopDes.begin(); desItr != inLoopDes.end(); ++desItr)
	{//for desItr
		inLoop.push_back(desItr->name);
	}//for desItr*/

	 // ueb forcing
	const char* ueb_forcNames[] = {
		"uebPrec",
		"uebTair",
		"uebTamin",
		"uebTamax",
		"uebVp",
		"uebWindS"
	};
	vector< string > inLoop;
	for (int ifn = 0; ifn < 6; ifn++)
	{
		inLoop.push_back(ueb_forcNames[ifn]);
	}
	//Ohd_Hydro::GridData1D::GridData1D(std::string outdir, std::vector< std::string > const& indIds, std::map< std::string, size_t > const& sibs,
	//	std::vector< std::string > const& vars, boost::posix_time::time_duration const& s, boost::posix_time::time_period prd, bool ignore1d);

	time_period_Ens = "20081001T07 20090830T07";

	GridData1D::GridData1D gd1d(any_cast<string> (deck.find("output-path")->second),
		get_property(g, independent_selected_basin_t()),
		getIndBasinSizes(g), 
		inLoop,
		get_property(g, graph_time_step_t()), //any_cast< time_duration >( deck[ "time-step" ] ),
		any_cast<time_period>(time_period_Ens),		//any_cast<time_period>(deck.find("time-period")->second),
		any_cast<bool>(deck.find("ignore-1d-xmrg")->second));

	//boost::posix_time::time_duration timeStep(get_property(allPixels, graph_time_step_t()));

	//boost::gregorian::date(year_type y, month_type m, day_type d)
	//ptime(gregorian::date d, time_duration_type td)
	boost::posix_time::ptime currentTime(boost::gregorian::date(year_in, month_in, day_in), (time_duration_type) td_in);				// (get_property(allPixels, graph_time_period_start_t()));

	//while (currentTime <= get_property(allPixels, graph_time_period_end_t())) {
	//get_property(allPixels, graph_time_t()) = currentTime;
	//readDataGrids< PropertyInitializer< InsideLoopPolicy > >(deck, allPixels, mptInsideLoop, gd1d);
	//pt::ptime curtime = get_property(allPixels, graph_time_t());

	typedef graph_property< PixelGraph, independent_selected_basin_t >::type IndSelBasinsType;
	IndSelBasinsType independentSelectedBasins = get_property(g, independent_selected_basin_t());
	IndSelBasinsType::iterator indItr;

	vector< float > val[6];
	for (int ifn = 0; ifn < 6; ifn++)
	{
		for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
		{//for indItr
			gd1d.readAtTime(currentTime, inLoop[ifn], *indItr, val[ifn]);			//gd1d.readAtTime(curtime, des.name, *indItr, val);
			//transform(val.begin(), val.end(), val.begin(), bind2nd(multiplies< float >(), std::abs(factor)));
			//int error = put_independent_basin_values(allPixels, indItr->c_str(), val, des);
		}//for indItr 
	}

	//currentTime = get_property(allPixels, graph_time_t()) +  get_property(allPixels, graph_time_step_t());  }

	cDateTime = cDateTime + get_property(g, graph_time_step_t());				// (get_property(allPixels, graph_time_period_start_t()));
	cout << " CurrentDatetime = " << cDateTime << endl;
	vector< float > val2[6];
	for (int ifn = 0; ifn < 6; ifn++)
	{
		for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
		{//for indItr
			gd1d.readAtTime(cDateTime, inLoop[ifn], *indItr, val2[ifn]);			//gd1d.readAtTime(curtime, des.name, *indItr, val);
			//transform(val.begin(), val.end(), val.begin(), bind2nd(multiplies< float >(), std::abs(factor)));
		    //int error = put_independent_basin_values(allPixels, indItr->c_str(), val, des);
		}//for indItr 
	}
	cout << " here -2 " << endl;
	for (int ifn = 0; ifn < 6; ifn++)
	{
		for (int ix = 0; ix < val2[ifn].size(); ix++)
			std::cout << " " << val2[ifn][ix];
		std::cout << std::endl;
		std::getchar();
	}
	std::getchar();
	cout << " here -8 " << endl;
	boost::posix_time::ptime cDateTimeT(boost::gregorian::date(CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2]),
		boost::posix_time::hours(8));
	cout << " CurrentDatetime = " << cDateTimeT << endl;
	vector< float > valT;
	for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
		gd1d.readAtTime(cDateTimeT, inLoop[1], *indItr, valT);
	for (int ix = 0; ix < valT.size(); ix++)
		std::cout << " " << valT[ix];
	std::cout << std::endl;
	cout << " here -9 " << endl;
	valT.clear();
	boost::posix_time::ptime cDateTimeT9(boost::gregorian::date(CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2]),
		boost::posix_time::hours(9));
	cout << " CurrentDatetime = " << cDateTimeT9 << endl;
	for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
		gd1d.readAtTime(cDateTimeT9, inLoop[1], *indItr, valT);
	for (int ix = 0; ix < valT.size(); ix++)
		std::cout << " " << valT[ix];
	std::cout << std::endl;
	cout << " here -10 " << endl;
	valT.clear();
	boost::posix_time::ptime cDateTimeT10(boost::gregorian::date(CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2]),
		boost::posix_time::hours(10));
	cout << " CurrentDatetime = " << cDateTimeT10 << endl;
	for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
		gd1d.readAtTime(cDateTimeT10, inLoop[1], *indItr, valT);
	for (int ix = 0; ix < valT.size(); ix++)
		std::cout << " " << valT[ix];
	std::cout << std::endl;





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

	// tamin and tamax one value per day --for the daArray from netCDF
	int currSimDay = (int)uebCellDAObjArr[0].istep / uebCellDAObjArr[0].nstepinaDay;
	if ((uebCellDAObjArr[0].istep % 120) == 0)
		std::cout << std::endl << " UEB time step: " << uebCellDAObjArr[0].istep << " current day: " << currSimDay << std::endl;

	//std::cout << " Obs R Matrix" << std::endl;	std::cout << objEnKFArr.R_obsErrCov << std::endl;	std::getchar();
	// get xmrg forcing "uebPrec",	"uebTair",	"uebTamin", "uebTamax", "uebVp", "uebWindS"
	float *uebForcArray = new float[6 * npix];
	if (error = get_values(g, ueb_forcNames, 6, uebForcArray)) return error;
	/*for (int ifn = 0; ifn < 6; ifn++)
	{
	for (int ix = 0; ix < npix; ix++)
	std::cout << " " << uebForcArray[ifn*npix + ix];
	std::cout << std::endl;
	std::getchar();
	}*/
	for (int ix = 0; ix < npix; ix++)
		std::cout << " " << uebForcArray[npix + ix];
	std::cout << std::endl;
	//std::cout << " Here -3: " << std::endl;
	int nOpts = objEnKFArr.numObsPoints;
	float* uebForcArrayDaPoint = new float[6 * nOpts];
	int retvalue = 0;
	uebCellDA *daUebCellObjArr = NULL;
	if (uebCellDAObjArr[0].daAssimlate)
	{
		daUebCellObjArr = new uebCellDA[npix];  // , objCell0);   // 10.6.17 cells (at obs points) for data assimilation
		if (error = get_value(g, "daUebCellObjArray", &daUebCellObjArr[0])) return error;
		//std::cout << " retrieved da_uebCellDAObjArray for time step " << uebCellDAObjArr[0].istep + 1 << std::endl;

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
		//std::cout << " Read forcing at obs points " <<  std::endl;
	}
	//std::cout << " Here -2: " << std::endl;
	//TODO---what if the assimilation obs. ends before simulation time end? 

	//start & end datetime
	int CurrentDateTime[4]; //check this	
	if (error = get_year_month_day_hour(g, CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2], CurrentDateTime[3]))	return error;
	double ModelDtMinutes, ModelDt;  //dt
	if (error = get_time_step_minutes(g, ModelDtMinutes)) return error;
	ModelDt = (double)ModelDtMinutes / 60.0;   //dT in hours
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

	 //cout << " here -10 " << endl;
	 // ueb forcing for ensemble run
	typedef graph_property< PixelGraph, independent_selected_basin_t >::type IndSelBasinsType;
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
	const char* time_period_Ens = "20081001T07,20090830T07";
	//cout << " here -7 " << endl;
	GridData1D gd1d(any_cast<string> (deck.find("output-path")->second),
		get_property(g, independent_selected_basin_t()),
		bSizes,				//getIndBasinSizes(g),		//std::map< std::string, size_t > const& sibs,
		inLoop,
		get_property(g, graph_time_step_t()), //any_cast< time_duration >( deck[ "time-step" ] ),
		any_cast<time_period>(deck.find("time-period")->second),		//any_cast<time_period>(time_period_Ens),		//
		any_cast<bool>(deck.find("ignore-1d-xmrg")->second));
	//cout << " here -6 " << endl;
	//cout << " ens period " << any_cast<time_period>(time_period_Ens)<<endl;
	//boost::posix_time::time_duration timeStep(get_property(allPixels, graph_time_step_t()));	
	boost::posix_time::ptime cDateTime(boost::gregorian::date(CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2]),
		boost::posix_time::hours(7)); // CurrentDateTime[3])); // //ptime(gregorian::date d, time_duration_type td)  //boost::gregorian::date(year_type y, month_type m, day_type d)
									  //boost::posix_time::ptime currentDateTime = currentDate +  get_property( g, graph_time_step_t() );				// (get_property(allPixels, graph_time_period_start_t()));
	cout << " CurrentDatetime = " << cDateTime << endl;
	//while (currentTime <= get_property(allPixels, graph_time_period_end_t())) {
	//get_property(allPixels, graph_time_t()) = currentTime;
	//readDataGrids< PropertyInitializer< InsideLoopPolicy > >(deck, allPixels, mptInsideLoop, gd1d);
	//pt::ptime curtime = get_property(allPixels, graph_time_t());
	//cout << " here -5 " << endl;
	vector< float > val[6];
	for (int ifn = 0; ifn < 6; ifn++)
	{
		for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
		{//for indItr
			gd1d.readAtTime(cDateTime, inLoop[ifn], *indItr, val[ifn]);			//gd1d.readAtTime(curtime, des.name, *indItr, val);
		}//for indItr 
	}
	/*cout << " here -4 " << endl;
	for (int ifn = 0; ifn < 6; ifn++)
	{
	for (int ix = 0; ix < val[ifn].size(); ix++)
	std::cout <<" "<<val[ifn][ix];
	std::cout << std::endl;
	}
	cout << " here -3 " << endl;
	for (int ix = 0; ix < val[1].size(); ix++)
	std::cout << " " << val[1][ix];
	std::cout << std::endl;*/

	vector< float > valT;
	InputDataType inData;
	if (deck.find("input-data") != deck.end())
	{//if
		inData = bt::any_cast<InputDataType>(deck.find("input-data")->second);
	}//if

	InputDataType::value_type::second_type::iterator pos;
	//InputDataType::const_iterator itr;
	InputDataType::iterator itri;

	ModelDataDescript tempDesc;
	if (error = get_data_descript(g, ueb_forcNames[1], tempDesc)) return error;  //get_data_descript(GraphType& g, const char* name, ModelDataDescript& des)

	for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
	{//for indItr

		gd1d.readAtTime(cDateTime, inLoop[1], *indItr, valT);

		float factor(-1.f);
		for (itri = inData.begin(); itri != inData.end(); ++itri)
		{//for itr
			if (itri->first == *indItr)
			{//if

				pos = itri->second.find(inLoop[1]);
				if (pos == itri->second.end())
				{

				}
				else
				{
					factor = pos->second;
				}
				break;
			}//if
		}//for itr
		if (itri == inData.end())
		{//if
			throw std::runtime_error("No such basin as " + *indItr + " in input deck!");
		}//if

		if (factor < 0.f && not_equal_to< float >()(factor, -1.f))
		{//if
			transform(valT.begin(), valT.end(), valT.begin(), bind2nd(multiplies< float >(), std::abs(factor)));
		}//if 

		error = put_independent_basin_values(g, indItr->c_str(), valT, tempDesc);

		if (error)
		{//if
			throw std::runtime_error("put_independent_basin_values Error!");
		}//if

		for (int ix = 0; ix < valT.size(); ix++)
			std::cout << " " << valT[ix];
		std::cout << std::endl;

	}//for indItr 

	float *uebT = new float[npix];
	if (error = get_value(g, ueb_forcNames[1], uebT)) return error;
	for (int ix = 0; ix < npix; ix++)
		std::cout << " " << uebT[ix];
	std::cout << std::endl;
	delete[] uebT;
	uebT = NULL;

	boost::posix_time::ptime cDateTime2(boost::gregorian::date(CurrentDateTime[0], CurrentDateTime[1], CurrentDateTime[2]),
		boost::posix_time::hours(10)); // CurrentDateTime[3])); // , (ptime::time_duration_type) CurrentDateTime[3]);				// (get_property(allPixels, graph_time_period_start_t()));
									   //boost::posix_time::ptime currentDateTime = currentDate +  get_property( g, graph_time_step_t() );				// (get_property(allPixels, graph_time_period_start_t()));
	cout << " CurrentDatetime = " << cDateTime2 << endl;
	vector< float > valT2;
	for (indItr = independentSelectedBasins.begin(); indItr != independentSelectedBasins.end(); ++indItr)
	{//for indItr

		gd1d.readAtTime(cDateTime2, inLoop[1], *indItr, valT2);

		float factor(-1.f);
		for (itri = inData.begin(); itri != inData.end(); ++itri)
		{//for itr
			if (itri->first == *indItr)
			{//if

				pos = itri->second.find(inLoop[1]);
				if (pos == itri->second.end())
				{

				}
				else
				{
					factor = pos->second;
				}
				break;
			}//if
		}//for itr
		if (factor < 0.f && not_equal_to< float >()(factor, -1.f))
		{//if
			transform(valT2.begin(), valT2.end(), valT2.begin(), bind2nd(multiplies< float >(), std::abs(factor)));
		}//if 

		error = put_independent_basin_values(g, indItr->c_str(), valT2, tempDesc);
		if (error)
		{//if
			throw std::runtime_error("put_independent_basin_values Error!");
		}//if

		for (int ix = 0; ix < valT2.size(); ix++)
			std::cout << " " << valT2[ix];
		std::cout << std::endl;
	}//for indItr 

	float *uebT2 = new float[npix];
	if (error = get_value(g, ueb_forcNames[1], uebT2)) return error;
	for (int ix = 0; ix < npix; ix++)
		std::cout << " " << uebT2[ix];
	std::cout << std::endl;
	delete[] uebT2;
	uebT2 = NULL;

	if (uebCellDAObjArr[0].istep == 2)
		exit(0);
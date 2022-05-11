#ifndef UEBDASRPF_H
#define UEBDASRPF_H

/*
* filename: uebDAsacrutpix7.h
* description: Define the UEB model with EnKF DA
*/

#include <iostream>
#include "ModelBase.hpp"
#include "PixelGraph.hpp"
#include "default_model_func.h"

using namespace Ohd_Hydro;

int uebDAsrPFFuncBeforeLoop(PixelGraph& g);
int uebDAsrPFFuncInsideLoop(PixelGraph& g);
//int uebDAsacrutpix7FuncAfterLoop(PixelGraph& g);

DEFINE_MODEL(uebDAsrPF,
	uebDAsrPFFuncBeforeLoop,
	uebDAsrPFFuncInsideLoop,
	default_model_func,  //
	default_model_func);

const ModelDataDescript uebDAsrPF::inputBeforeLoop[] = {

//11.29.17 UEB parameters to be read from deck as user-input

//*********UEB Site Variables---distributed params*******************
	{
		"ueb_df", //Drift factor multiplier
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_apr", //Average atmospheric pressure 
		"Pa",        //unit
		"FL^-2",         //dimension
		"MAPX",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_Aep", // Albedo extinction coefficient 
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_cc", // Canopy coverage fraction   
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_hcan", // Canopy height  
		"m",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_lai", // Leaf area index
		"m^2/m^2",        //unit
		"L^2/L^2",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_Sbar", // Maximum snow load held per unit branch area 
		"kg/m^2",        //unit
		"M/L^2",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_ycage", // Forest age flag for wind speed profile parameterization  
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0,                 //missing value
		0,                 //missing value for the first time step
		-1,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_slope", // A 2-D grid that contains the slope at each grid point  
		"Deg",        //unit
		"Deg",         //dimension
		"MAPX",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_aspect", // A 2-D grid that contains the aspect at each grid point 
		"Deg",        //unit
		"Deg",         //dimension
		"MAPX",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	
	{
		"ueb_subalb", // Albedo (fraction 0-1) of the substrate beneath the snow (ground, or glacier)
		"N/A",        //unit
			"DLESS",         //dimension
			"MAPX",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_subtype", // Type of beneath snow substrate encoded as (0 = Ground/Non Glacier, 1=Clean Ice/glacier, 2= Debris covered ice/glacier, 3= Glacier snow accumulation zone)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0,                 //missing value
			0,                 //missing value for the first time step
			-1,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_gsurf", // The fraction of surface melt that runs off (e.g. from a glacier)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b01", // Bristow-Campbell B for January (1)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b02", // Bristow-Campbell B for February (2)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b03", // Bristow-Campbell B for March(3)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b04", // Bristow-Campbell B for April (4)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b05", // Bristow-Campbell B for may (5)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b06", // Bristow-Campbell B for June (6)
		"N/A",        //unit
			"DLESS",         //dimension
			"SQIN",         //Data Type Code
			UNKNOWN,     //accumulation policy
			0,           //max number of missing
			DONT_FILLMISSING,    //if fill missing
			0.f,                 //missing value
			0.f,                 //missing value for the first time step
			-1.f,                 // no data value
			NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b07", // Bristow-Campbell B for July (7)
		"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data Type Code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b08", // Bristow-Campbell B for August (8)
		"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data Type Code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b09", // Bristow-Campbell B for September (9)
		"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data Type Code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b10", // Bristow-Campbell B for October (10) 
		"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data Type Code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b11", // Bristow-Campbell B for November (11)
		"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data Type Code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
	},
	{
		"ueb_b12", // Bristow-Campbell B for December (12)
		"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data Type Code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
	},
	//
	// Qg as site var 11.27.17
	{
		"ueb_Qg", //ground energy flux
		"KJ/m^2-hr",        //unit
		"M/L^2-T",         //dimension
		"SQIN",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
	// SnowAlb as site var 11.27.17
	{
		"ueb_SnowAlb", //ground energy flux 
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
    //
//*********UEB states************************************************************ 
	//
	{
		"Us", //Energy content
		"KJ/m^2",        //unit
		"M/L^2",         //dimension
		"SQIN",         //Data Type Code
		SUM,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
	{
		"Ws", //Snow water equivalent
		"m",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
	{
		"Wc", //Snow water equivalent canopy
		"m",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
	{
		"Tausn", //Snow surface dimensionless age
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		0.f,                 //missing value
		0.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},

//SAC input before 
//
//
// SAC_SMA parameters
//
//
// first SAC_SMA parameter
//
	{
		"sac_UZTWM", //name
		"MM",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// second SAC_SMA parameter
		//
	{
		"sac_UZFWM", //name
		"MM",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 3rd SAC_SMA parameter
		//
	{
		"sac_UZK", //name
		"1/deg",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 4th SAC_SMA parameter
		//
	{
		"sac_PCTIM", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 5th SAC_SMA parameter
		//
	{
		"sac_ADIMP", //name
		"MM",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 6th SAC_SMA parameter
		//
	{
		"sac_RIVA", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 7th SAC_SMA parameter
		//
	{
		"sac_ZPERC", //name
		"MM",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 8th SAC_SMA parameter
		//
	{
		"sac_REXP", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 9th SAC_SMA parameter
		//
	{
		"sac_LZTWM", //name
		"MM",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 10th SAC_SMA parameter
		//
	{
		"sac_LZFSM", //name
		"MM",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 11th SAC_SMA parameter
		//
	{
		"sac_LZFPM", //name
		"MM",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 12th SAC_SMA parameter
		//
	{
		"sac_LZSK", //name
		"1/deg",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 13th SAC_SMA parameter
		//
	{
		"sac_LZPK", //name
		"1/deg",        //unit
		"L",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 14th SAC_SMA parameter
		//
	{
		"sac_PFREE", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 15th SAC_SMA parameter
		//
	{
		"sac_SIDE", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 16th SAC_SMA parameter
		//
	{
		"sac_RSERV", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// 17th SAC_SMA parameter --- fraction of forest cover
		//
	{
		"sac_EFC", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		//PE
		//
		//
		// Jan PE
		//
	{
		"pe_JAN", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// Feb PE
		//
	{
		"pe_FEB", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// MAR PE
		//
	{
		"pe_MAR", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// APR PE
		//
	{
		"pe_APR", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// MAY PE
		//
	{
		"pe_MAY", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// JUN PE
		//
	{
		"pe_JUN", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// JUL PE
		//
	{
		"pe_JUL", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// AUG PE
		//
	{
		"pe_AUG", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// SEP PE
		//
	{
		"pe_SEP", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// OCT PE
		//
	{
		"pe_OCT", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// NOV PE
		//
	{
		"pe_NOV", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// DEC PE
		//
	{
		"pe_DEC", //name
		"mm/day",        //unit
		"L/T",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		//PE Adj
		//
		//
		// Jan PE adj
		//
	{
		"peadj_JAN", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// Feb PE adj
		//
	{
		"peadj_FEB", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// MAR PE adj
		//
	{
		"peadj_MAR", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// APR PE adj
		//
	{
		"peadj_APR", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// MAY PE adj
		//
	{
		"peadj_MAY", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// JUN PE adj
		//
	{
		"peadj_JUN", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// JUL PE adj
		//
	{
		"peadj_JUL", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// AUG PE adj
		//
	{
		"peadj_AUG", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// SEP PE adj
		//
	{
		"peadj_SEP", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// OCT PE adj
		//
	{
		"peadj_OCT", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// NOV PE adj
		//
	{
		"peadj_NOV", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// DEC PE adj
		//
	{
		"peadj_DEC", //name
		"N/A",        //unit
		"DLESS",         //dimension
		"SQIN",         //Data Type Code
		UNKNOWN,     //accumulation policy
		0,           //max number of missing
		DO_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		NO_TIME_INTERVAL     //time interval
	},
		//
		// SAC_SMA states in precentage
		//
		//
		// 1st SAC_SMA state
		//
	{
		"uztwc", //name
		"PCTD",        //unit
		"DLESS",         //dimension
		"SASC",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 2nd SAC_SMA state
		//
	{
		"uzfwc", //name
		"PCTD",        //unit
		"DLESS",         //dimension
		"SASC",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 3rd SAC_SMA state
		//
	{
		"lztwc", //name
		"PCTD",        //unit
		"DLESS",         //dimension
		"SASC",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 4th SAC_SMA state
		//
	{
		"lzfsc", //name
		"PCTD",        //unit
		"DLESS",         //dimension
		"SASC",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 5th SAC_SMA state
		//
	{
		"lzfpc", //name
		"PCTD",        //unit
		"DLESS",         //dimension
		"SASC",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 6th SAC_SMA state
		//
	{
		"adimpc", //name
		"PCTD",        //unit
		"DLESS",         //dimension
		"SASC",         //Data Type Code
		AVERAGE,     //accumulation policy
		0,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},

		//
		// SAC_SMA states in real value not precentage
		//
		//
		// 1st SAC_SMA state
		//
	{
		"real_uztwc", //name
		"MM",        //unit
		"L",         //dimension
		"MAP",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 2nd SAC_SMA state
		//
	{
		"real_uzfwc", //name
		"MM",        //unit
		"L",         //dimension
		"MAP",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 3rd SAC_SMA state
		//
	{
		"real_lztwc", //name
		"MM",        //unit
		"L",         //dimension
		"MAP",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 4th SAC_SMA state
		//
	{
		"real_lzfsc", //name
		"MM",        //unit
		"L",         //dimension
		"MAP",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 5th SAC_SMA state
		//
	{
		"real_lzfpc", //name
		"MM",        //unit
		"L",         //dimension
		"MAP",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
		//
		// 6th SAC_SMA state
		//
	{
		"real_adimpc", //name
		"MM",        //unit
		"L",         //dimension
		"MAP",         //Data Type Code
		AVERAGE,     //accumulation policy
		UNLIMITED_MISSING,   //max number of missing
		DONT_FILLMISSING,    //if fill missing
		-1.f,                 //missing value
		-1.f,                 //missing value for the first time step
		-1.f,                 // no data value
		"1:00:00"     //time interval
	},
//
//
// rutpix7 parameters
//
//
		// first rutpix7  parameter
		//
		{
			"rutpix_SLOPC", //name
				"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data type code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
			//
			// second rutpix7 parameter
			//
		{
			"rutpix_ROUGC", //name
			"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data type code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 3rd rutpix7 parameter, input as mm/hr but 'do_snow' converts to 
				//                       mm/dthr
				//
		{
			"rutpix_BETAC", //name
			"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data type code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 4th rutpix7 parameter, input as mm/hr but 'do_snow' converts to 
				//                       mm/dthr
				//
		{
			"rutpix_ALPHC", //name
			"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data type code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 5th rutpix7 parameter
				//
		{
			"rutpix_SLOPH", //name
			"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data type code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 6th rutpix7 parameter
				//
		{
			"rutpix_DS", //name
			"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data type code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 7th rutpix7 parameter
				//
		{
			"rutpix_ROUGH", //name
			"N/A",        //unit
				"DLESS",         //dimension
				"SQIN",         //Data type code
				UNKNOWN,     //accumulation policy
				0,           //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 8th rutpix7 parameter
				//
		{
			"rutpix_qlos_cnst",           //name
			"CMS",                        //unit
				"M3/T",                       //dimension
				"SQIN",                       //Data type code
				UNKNOWN,                      //accumulation policy
				UNLIMITED_MISSING,            //max number of missing
				DONT_FILLMISSING,             //if fill missing
				0.f,                          //missing value
				0.f,                          //missing value for the first time step
				-1.f,                          // no data value
				NO_TIME_INTERVAL               //time interval
		},
				//
				// 9th rutpix7 parameter
				//
		{
			"rutpix_qlos_rate",           //name
			"N/A",                        //unit
				"DLESS",                       //dimension
				"SQIN",                       //Data type code
				UNKNOWN,                      //accumulation policy
				UNLIMITED_MISSING,            //max number of missing
				DONT_FILLMISSING,             //if fill missing
				-999.f,                       //missing value
				-999.f,                       //missing value for the first time step
				-1.f,                          // no data value
				NO_TIME_INTERVAL               //time interval
		},
				//
				// 10th rutpix7 parameter
				//
		{
			"rutpix_qlos_pow",           //name
			"N/A",                        //unit
				"DLESS",                       //dimension
				"SQIN",                       //Data type code
				UNKNOWN,                      //accumulation policy
				UNLIMITED_MISSING,            //max number of missing
				DONT_FILLMISSING,             //if fill missing
				-999.f,                       //missing value
				-999.f,                       //missing value for the first time step
				-1.f,                          // no data value
				NO_TIME_INTERVAL               //time interval
		},
				//
				// rutpix7 states 
				//
				//
				// 1st rutpix7 state
				//
		{
			"areac", //name
			"m^2",        //unit
				"L^2",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				"1:00:00"     //time interval
		},
				//
				// 2nd rutpix7 state
				//
		{
			"areac1", //name
			"m^2",        //unit
				"L^2",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				"1:00:00"     //time interval
		},
				//
				// 3rd rutpix7 state
				//
		{
			"areac2", //name
			"m^2",        //unit
				"L^2",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				"1:00:00"     //time interval
		},
				//
				// 4th rutpix7 state
				//
		{
			"areac3", //name
			"m^2",        //unit
				"L^2",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				"1:00:00"     //time interval
		},
				//
				// 5th rutpix7 state
				//
		{
			"areac4", //name
			"m^2",        //unit
				"L^2",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				"1:00:00"     //time interval
		},
				//
				// 6th rutpix7 state
				//
		{
			"depth", //name
			"mm",        //unit
				"L",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				"1:00:00"     //time interval
		},
				//
				// 1st rutpix7 uhg
				//
		{
			"uhg1", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 2nd rutpix7 uhg
				//
		{
			"uhg2", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 3rd rutpix7 uhg
				//
		{
			"uhg3", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 4th rutpix7 uhg
				//
		{
			"uhg4", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 5th rutpix7 uhg
				//
		{
			"uhg5", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 6th rutpix7 uhg
				//
		{
			"uhg6", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 7th rutpix7 uhg
				//
		{
			"uhg7", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 8th rutpix7 uhg
				//
		{
			"uhg8", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 9th rutpix7 uhg
				//
		{
			"uhg9", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 10th rutpix7 uhg
				//
		{
			"uhg10", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				-1.f,                 //missing value
				-1.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// 11th rutpix7 uhg
				//
		{
			"nord", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				NO_TIME_INTERVAL     //time interval
		},
				//
				// channel loss state
				//
		{
			"channelloss", //name
			"N/A",        //unit
				"N/A",         //dimension
				"SQIN",         //Data type code
				AVERAGE,     //accumulation policy
				UNLIMITED_MISSING,   //max number of missing
				DONT_FILLMISSING,    //if fill missing
				0.f,                 //missing value
				0.f,                 //missing value for the first time step
				-1.f,                 // no data value
				"1:00:00"     //time interval
		}



};
const int uebDAsrPF::numOfInputBeforeLoop = 112;   // 31 + 53 + 28;

const ModelDataDescript uebDAsrPF::inputInsideLoop[] = {
	{
		"uebPrec", 
		"m/hr", 
	    "L/T", 
	    "uebPrec",
	    AVERAGE,
		0,
		DONT_FILLMISSING,
		0.f,
	    0.f,
		-9999.f,
		"1:00:00"
	},
	{
		"uebTair", 
		"DEGC", 
	    "T", 
	    "uebTair",
		AVERAGE,
		0,
		DONT_FILLMISSING,
		0.f,
		0.f,
		-9999.f,
		"1:00:00"
	},
	{
		"uebTamin", 
		"DEGC",
		"T",
		"uebTamin",
		AVERAGE,
		23,
		DO_FILLMISSING,
	    USE_PREVIOUS,
	    -4.f,
		-9999.f,
		"24:00:00"
	},
	{
		"uebTamax",
		"DEGC",
		"T",
		"uebTamax",
		AVERAGE,
		23,
		DO_FILLMISSING,
		USE_PREVIOUS,
		4.f,
		-9999.f,
		"24:00:00"
	},
	{
		"uebVp",
		"Pa", 
		"F/L^2",  
		"uebVp",
		AVERAGE,
		0,
		DONT_FILLMISSING,
		400.f,
		400.f,
		-9999.f,
		"1:00:00"
	},
	{
		"uebWindS", 
		"m/s",
		"L/T",  
		"uebWindS",
		AVERAGE,
		0,
		DONT_FILLMISSING,
		2.f,
		2.f,
		-9999.f,
		"1:00:00"
	},

	{
		"modisDelVis",
		"%",
		"%",
		"modisDelVis",
		AVERAGE,
		UNLIMITED_MISSING,
		DONT_FILLMISSING,
		USE_PREVIOUS,
		0.f,
		-1.f,
		"24:00:00"
	},

//7.13.18 using RDHM-SACSMA units for Temp and prec (F, and mm / dT)
/*	{
		"uebPrec", "MM", "L", "MAPX",                  //xmrg
		SUM,
		UNLIMITED_MISSING,
		DONT_FILLMISSING,
		0.f,
		0.f,
		-1.f,
		"1:00:00"
	},
	{
		"uebTair", "DEGF", "T", "TEMP",				//tair
		AVERAGE,
		12,
		DO_FILLMISSING,
		40.f,
		USE_PREVIOUS,
		-99.f,
		"1:00:00"
	},
	{
		"uebTamin", "DEGF", "T", "TEMP",
		AVERAGE,
		24,
		DO_FILLMISSING,
		USE_PREVIOUS,
		24.f,
		-99.f,
		"24:00:00"
	},
	{
		"uebTamax", "DEGF", "T", "TEMP",
		AVERAGE,
		24,
		DO_FILLMISSING,
		USE_PREVIOUS,
		40.f,
		-99.f,
		"24:00:00"
	},
*/
// No SAC inside-loop input--xmrg coming from UEB

//                                 	"tair", "C", "T",
//	                                 AVERAGE,
//						       12,
//					 DO_FILLMISSING,
//						       40.f,
//					 USE_PREVIOUS,
//						       -99.f,
//						       "1:00:00"
	/*{ "xmrg", "MM", "L", "MAPX",
	SUM,
	UNLIMITED_MISSING,
	DONT_FILLMISSING,
	0.f,
	0.f,
	-1.f,
	"1:00:00" },*/
	/*       {
	"surf_water", //name
	"N/A",        //unit
	"DLESS",         //dimension
	"MAP",         //Data Type Code
	AVERAGE,     //accumulation policy
	UNLIMITED_MISSING,   //max number of missing
	DONT_FILLMISSING,    //if fill missing
	USE_PREVIOUS,                 //missing value
	-1.f,                 //missing value for the first time step
	-1.f,                 // no data value
	"1:00:00"     //time interval
	}
	*/

// No rutpix7 inside input

};

const int uebDAsrPF::numOfInputInsideLoop = 7;      //  = 6 + 0 + 0

const ModelDataDescript uebDAsrPF::inputAfterLoop[] = {};

const int uebDAsrPF::numOfInputAfterLoop = 0;       //0 + 0 + 0

#endif//#ifndef UEBDASRPF_H


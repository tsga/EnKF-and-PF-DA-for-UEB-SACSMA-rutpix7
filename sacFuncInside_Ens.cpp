/*
* filename : sacFuncInside_Ens.cpp

4.25.18
A stripped down SAC call (with no frz and snow) for ensemble simulation of Q with UEB

*/
#include <vector>
#include <iostream>
#include <iterator>
#include <boost/property_map/property_map.hpp>
#include <boost/shared_array.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "hlrms_interfaces.hpp"
#include "getFromDeck.hpp"
#include "getSelectedVerticesInConnOrder.hpp"

#include "uebDAsrPFuebpgdecls.h"

extern "C" {
#include "models.h"
#ifdef CHPS
#include "utilities.h"
#include "logging.h"
#include "fortranCommonIncludes.h"
#endif
}//extern "C"

extern const Ohd_Hydro::DECK deck;

using namespace Ohd_Hydro;

using namespace std;

//
//void do_sac_EnsHost_PC(int npix, SACGridDA *sacGridArray);
//int do_sac_EnsHost(int threadsPerBlock, int npix, SACGridDA *sacGridArray);

int sacFuncInside_Ens(PixelGraph& g, int threadsPerBlock, std::vector<float> pcp, SACGridDA *sacGridArray, int gpu0)
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

	int npix;
	int error = retrieve_value(g, "npix", npix);
	std::vector<int> do_sac_error(npix, 0);

	typedef boost::graph_traits< PixelGraph >::vertex_descriptor Vertex;
	bt::shared_array< Vertex > selectedPixels;
	error = retrieve_value(g, "selectedPixels", selectedPixels);
	//cout << " here in_1 " << endl;
	typedef boost::property_map < PixelGraph, attribute_t >::type AttMapType;
	AttMapType  attMap = boost::get(attribute_t(), g);
	typedef boost::graph_property< PixelGraph, name_index_t >::type NameIndexType;
	NameIndexType& nameIndex = boost::get_property(g, name_index_t());
	NameIndexType::iterator namePos;

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
	//cout << " here in_2 " << endl;
	for (int ix = 0; ix < npix; ix++)
	{
		sacGridArray[ix].dtm = (float)dtm;
		sacGridArray[ix].dtday = (float)dtm / (24. * 60);  //convert time step in days
		sacGridArray[ix].aesc = 0;
		for (size_t i = 0; i < 17; ++i)
			sacGridArray[ix].sacpar[i] = sac_par[ix * 17 + i];
	}

	ModelDataDescript pcpDes;
	if (boost::get_property(g, data_descript_t()).find("xmrg") != boost::get_property(g, data_descript_t()).end())
	{//if
		pcpDes = boost::get_property(g, data_descript_t())["xmrg"];
	}//if
	else
	{//else
		std::cerr << "xmrg has not been defined in the descritor" << std::endl;
		return  PROPERTY_NOT_FOUND;
	}//else


	/* VK 05/17/07 Introduced averaging time interval, dt_avg  */
	static int count = 0;
	int dt_avg = 24;
	/*       int dt_avg = 1;    */
	if (count % (dt_avg * 60 / dtm) == 0)
	{
		count = 0;
	}
	//cout << " here in_5 " << endl;

	HRAP< float > loc;
	size_t pix(0);
	//
	for (Vertex *itr = selectedPixels.get(); itr < selectedPixels.get() + npix; ++itr)
	{//for itr
		namePos = nameIndex.find("Location");
		loc = boost::any_cast<HRAP<float>>(attMap[*itr][namePos->second]);
		sacGridArray[pix].hrapx = loc.x;
		sacGridArray[pix].hrapy = loc.y;

		get_ped(year, month, day, 1, dtm, pe.get() + 12 * pix, pe_adj.get() + 12 * pix, &sacGridArray[pix].edmnd);

		namePos = nameIndex.find("ped");
		if (namePos == nameIndex.end())
		{//if
			std::cerr << "Propert: ped" << " not found, Pixel =" << *itr << std::endl;
			return PROPERTY_NOT_FOUND;
		}//if
		attMap[*itr][namePos->second] = boost::any(sacGridArray[pix].edmnd);


		for (size_t i = 0; i < 6; ++i)
		{//for i
			namePos = nameIndex.find(sac_real_stnames[i]);
			if (namePos == nameIndex.end())
			{//if
				std::cerr << "Propert: " << sac_real_stnames[i] << " not found, Pixel =" << *itr << std::endl;
				return PROPERTY_NOT_FOUND;
			}//if
			 //sac_st_real.push_back(boost::any_cast<float>(attMap[*itr][namePos->second]));
			sacGridArray[pix].sacst[i] = boost::any_cast<float>(attMap[*itr][namePos->second]);

		}//for i
		nameIndex.find("xmrg");
		if (namePos == nameIndex.end())
		{//if

			std::cerr << "Property: xmrg/rmlt not found" << std::endl;
			return PROPERTY_NOT_FOUND;
		}//if
		 //float pcp = boost::any_cast< float >( attMap[ *itr ][ namePos->second ] );
		sacGridArray[pix].pxv = pcp[pix];

		if (pcp[pix] < pcpDes.noDataValue)   //if ( pcp[pix] <= pcpDes.noDataValue ) 
		{//if
			pcp[pix] = pcpDes.missingValue;
			attMap[*itr][namePos->second] = boost::any(pcp[pix]);
			sacGridArray[pix].pxv = pcp[pix];
		}//if

		//tet is output
		namePos = nameIndex.find("tet");
		if (namePos != nameIndex.end())
		{//if
			sacGridArray[pix].tet = boost::any_cast<float>(attMap[*itr][namePos->second]);
		}//if
		else sacGridArray[pix].tet = 0.f;

		if (sacGridArray[pix].pxv < 0.0)
		{
			sacGridArray[pix].pxv = 0.0;
			cout << "WARNING: Negative precip found.  Setting to zero.";
		}
		//cout << " here in_5 " << endl;	

		if (count == 0)
		{
			sacGridArray[pix].surf = 0.f;
			sacGridArray[pix].grnd = 0.f;
			sacGridArray[pix].tet = 0.f;
		}

		pix++;

	} //for itr

	//std::cout << " pix = " << pix << "count = " << count << endl;
	if (gpu0 == 0)
	{
		//cout << " here in_6a size of sacst_prv = " << sacst_prv.size();
		do_sac_EnsHost(threadsPerBlock, npix, sacGridArray);
	}
	else
	{
		do_sac_EnsHost_PC(npix, sacGridArray);
	}

	//
	pix = 0;
	//
	for (Vertex *itr = selectedPixels.get(); itr < selectedPixels.get() + npix; ++itr)
	{//for itr

		if (sacGridArray[pix].error == 1)
		{
			std::cout << "ERROR: operation: ""sac"" !!!! " << std::endl;

			std::cout << " at Pixel = ( " << sacGridArray[pix].hrapx << ", " << sacGridArray[pix].hrapy << " )" 
				<< "at step (YYYY/MM/DD/HH): " << year << '/' << month << '/' << day << '/' << hour << std::endl;

			std::cout << "Precipitation = " << pcp[pix] << " MM" << endl;
			
			/*std::cout << "SAC PAR = " << endl;
			copy(sac_par.get() + pix * 17, sac_par.get() + pix * 17 + 17, std::ostream_iterator< float >(std::cout, " "));
			std::cout << endl;
			std::cout << "PE = " << endl;
			copy(pe.get() + pix * 12, pe.get() + pix * 12 + 12, ostream_iterator< float >(std::cout, " "));
			std::cout << endl;
			std::cout << "PE_adj = " << endl;
			copy(pe_adj.get() + pix * 12, pe_adj.get() + pix * 12 + 12, ostream_iterator< float >(std::cout, " "));*/

			std::cout << "SAC ST = " << endl;
			for (int ist = 0; ist < 6; ist++)
				std::cout << "  " << sacGridArray[pix].sacst[ist];
			std::cout << endl;

			return ERROR;
		}

		vector< float > sac_st(6);
		float	  totalWater1 = 0.f;
		float	  totalWater2 = 0.f;

		sacGridArray[pix].updateSACGrid(sac_st, totalWater1, totalWater2);

		//update state values
		for (size_t i = 0; i < 6; ++i)
		{//for i
			attMap[*itr][nameIndex[sac_real_stnames[i]]] = boost::any(sacGridArray[pix].sacst[i]);
			attMap[*itr][nameIndex[sac_stnames[i]]] = boost::any(sac_st[i]);
		}//for i

		 //#ifdef SUB_SURF_FLOW
		attMap[*itr][nameIndex["surfaceFlow"]] = boost::any(sacGridArray[pix].surf);
		attMap[*itr][nameIndex["subsurfaceFlow"]] = boost::any(sacGridArray[pix].grnd);
		//#endif //#ifdef SUB_SURF_FLOW		
		attMap[*itr][nameIndex["tet"]] = boost::any(sacGridArray[pix].tet);

		attMap[*itr][nameIndex["totalWater1"]] = boost::any(totalWater1);
		attMap[*itr][nameIndex["totalWater2"]] = boost::any(totalWater2);

		attMap[*itr][nameIndex["interflow"]] = boost::any(sacGridArray[pix].sif);
		attMap[*itr][nameIndex["primaryFlow"]] = boost::any(sacGridArray[pix].bfp);
		attMap[*itr][nameIndex["supplementalFlow"]] = boost::any(sacGridArray[pix].bfs);
		attMap[*itr][nameIndex["excessFlow"]] = boost::any(sacGridArray[pix].ssur);
		attMap[*itr][nameIndex["directFlow"]] = boost::any(sacGridArray[pix].sdro);
		attMap[*itr][nameIndex["percolation"]] = boost::any(sacGridArray[pix].sperc);

		pix++;

	}//for itr

	count++;

	return OK;
}//sacFuncInside_Ens

__host__ __device__
void SACGridDA::updateSACGrid
(std::vector< float > &sac_st, float &totalWater1, float &totalWater2)
{
	//convert real state to precentage 
	//sac_st.resize(6);
	sac_st[0] = sacst[0] / sacpar[0];
	sac_st[1] = sacst[1] / sacpar[1];
	sac_st[2] = sacst[2] / sacpar[8];
	sac_st[3] = sacst[3] / sacpar[9];
	sac_st[4] = sacst[4] / sacpar[10];
	sac_st[5] = sacst[5] / (sacpar[0] + sacpar[8]);

	totalWater1 = (sacst[0] + sacst[1]) / (sacpar[0] + sacpar[1]);
	totalWater2 = (sacst[2] + sacst[3] + sacst[4]) / (sacpar[8] + sacpar[9] + sacpar[10]);
	
	return;
}

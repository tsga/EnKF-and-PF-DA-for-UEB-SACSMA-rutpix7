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

//--without-mpi
//2.25.18 writes vector to a 2D nc at time step
__device__ __host__
int Write1DVector_to2DNC(const char* FileName, const char* varName, int tIndx, int Nz_dim, float* var_inp) //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, v_varid = 0;

	const int  NDIMS = 2;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Open netcdf file.  
	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))             //| NC_MPIIO, inpComm, inpInfo
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, varName, &v_varid)))
		ERR(retncval);

	start[0] = tIndx;	//0;
	start[1] = 0;   // z_dim
	count[0] = 1;   // Nt_dim;
	count[1] = Nz_dim;

	//put variable values
	if (retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[0]))
		ERR(retncval);

	//close files
	if (retncval = nc_close(ncid))
		ERR(retncval);

	//std::cout << "Sucess storing dimension vars in: " << FileName << std::endl;
	//fflush(stdout); 
	return 0;
}
//--without-mpi
//1.25.18 writes a 2D array at time step
__device__ __host__
int Write2DSlub_to3DNC(const char* FileName, const char* varName, int tIndx, int Np_dim, int Nz_dim, float** var_inp) //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, v_varid = 0;

	const int  NDIMS = 3;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Open netcdf file.  
	if ((retncval = nc_open(FileName, NC_WRITE, &ncid)))             //| NC_MPIIO, inpComm, inpInfo
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, varName, &v_varid)))
		ERR(retncval);

	start[0] = tIndx;	//0;
	start[1] = 0;   // point
	start[2] = 0;   // z_dim;
	count[0] = 1;   // Nt_dim;
	count[1] = Np_dim;
	count[2] = Nz_dim;
		
	//put variable values
	if (retncval = nc_put_vara_float(ncid, v_varid, start, count, &var_inp[0][0]))
		ERR(retncval);
		
	//close files
	if (retncval = nc_close(ncid))
		ERR(retncval);

	//std::cout << "Sucess storing dimension vars in: " << FileName << std::endl;
	//fflush(stdout); 
	return 0;
}
//no-mpi: 2.25.18 for ens and da outputs 
//creates 2D netcdf and stores dimension variables; called once for a given output netcdf 
__device__ __host__
int Create2DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
	const char* zName, int tDim, int zDim, float* t_inp, float* fillVal)           // MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, t_dimid = 0, z_dimid = 0;
	int v_varid = 0, t_varid = 0;
	const char* vunits = "UNITS";
	const char* fillValue = "_FillValue";
	//float missVal = -9999;
	const int  NDIMS = 2;
	int dimids[NDIMS];
	int oldFill = 0;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Error handling.  
	int retncval = 0;
	// Create netcdf file.  
	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid)))
		ERR(retncval);
	//set fill on
	if ((retncval = nc_set_fill(ncid, NC_FILL, &oldFill)))
		ERR(retncval);
	/* Define the dimensions. record dim can be unlimited*/
	if ((retncval = nc_def_dim(ncid, tName, tDim, &t_dimid)))
		ERR(retncval);
	if ((retncval = nc_def_dim(ncid, zName, zDim, &z_dimid)))
		ERR(retncval);
	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
	dimids[0] = t_dimid;
	dimids[1] = z_dimid;
	// Define the netCDF variables 
	if ((retncval = nc_def_var(ncid, VarName, NC_FLOAT, NDIMS, dimids, &v_varid)))
		ERR(retncval);
	//assign missing
	if ((retncval = nc_put_att_float(ncid, v_varid, fillValue, NC_FLOAT, 1, fillVal)))
		ERR(retncval);
	// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, v_varid, vunits, strlen(varUnits), varUnits)))
		ERR(retncval);
	if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1, &t_dimid, &t_varid)))
		ERR(retncval);
	// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, t_varid, vunits, strlen(tUnitsout), tUnitsout)))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid, t_varid, "long_name", strlen(tlong_name), tlong_name))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid, t_varid, "calendar", strlen(tcalendar), tcalendar))
		ERR(retncval);

	//put values to dim variables
	if ((retncval = nc_put_var_float(ncid, t_varid, &t_inp[0])))
		ERR(retncval);

	//close file
	if ((retncval = nc_close(ncid)))
		ERR(retncval);
	//delte 3D array	
	std::cout << "Sucess creating and storing dimension vars in: " << FileName << std::endl;
	//fflush(stdout); 
	return 0;
}
//no-mpi: 1.18.18 for ens and da outputs 
//creates 3D netcdf and stores dimension variables; called once for a given output netcdf 
__device__ __host__
int Create3DNC_uebOutputs(const char* FileName, const char* VarName, const char *varUnits, const char* tName, const char* tUnitsout, const char* tlong_name, const char* tcalendar,
	                                                       const char* yxName, const char* zName, int tDim, int yxDim, int zDim, float* t_inp, float* fillVal)           // MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	// IDs for the netCDF file, dimensions, and variables. 
	int ncid = 0, t_dimid = 0, yx_dimid = 0, z_dimid = 0;
	int v_varid = 0, t_varid = 0;
	const char* vunits = "UNITS";
	const char* fillValue = "_FillValue";
	//float missVal = -9999;
	const int  NDIMS = 3;
	int dimids[NDIMS];
	int oldFill = 0;
	// The start and count arrays will tell the netCDF library where to write our data.    
	size_t start[NDIMS], count[NDIMS];
	// Error handling.  
	int retncval = 0;
	// Create netcdf file.  
	if ((retncval = nc_create(FileName, NC_NETCDF4 | NC_CLOBBER, &ncid)))
		ERR(retncval);
	//set fill on
	if ((retncval = nc_set_fill(ncid, NC_FILL, &oldFill)))
		ERR(retncval);
	/* Define the dimensions. record dim can be unlimited*/
	if ((retncval = nc_def_dim(ncid, tName, tDim, &t_dimid)))
		ERR(retncval);
	if ((retncval = nc_def_dim(ncid, yxName, yxDim, &yx_dimid)))
		ERR(retncval);
	if ((retncval = nc_def_dim(ncid, zName, zDim, &z_dimid)))
		ERR(retncval);
	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. In C, the unlimited dimension must come first on the list of dimids. */
	dimids[0] = t_dimid;
	dimids[1] = yx_dimid;
	dimids[2] = z_dimid;
	// Define the netCDF variables 
	if ((retncval = nc_def_var(ncid, VarName, NC_FLOAT, NDIMS, dimids, &v_varid)))
		ERR(retncval);
	//assign missing
	if ((retncval = nc_put_att_float(ncid, v_varid, fillValue, NC_FLOAT, 1, fillVal)))
		ERR(retncval);
	// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, v_varid, vunits, strlen(varUnits), varUnits)))
		ERR(retncval);
	if ((retncval = nc_def_var(ncid, tName, NC_FLOAT, 1, &t_dimid, &t_varid)))
		ERR(retncval);
	// Assign units attributes to the netCDF variables.  
	if ((retncval = nc_put_att_text(ncid, t_varid, vunits, strlen(tUnitsout), tUnitsout)))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid, t_varid, "long_name", strlen(tlong_name), tlong_name))
		ERR(retncval);
	if (retncval = nc_put_att_text(ncid, t_varid, "calendar", strlen(tcalendar), tcalendar))
		ERR(retncval);

	//put values to dim variables
	if ((retncval = nc_put_var_float(ncid, t_varid, &t_inp[0])))
		ERR(retncval);

	//close file
	if ((retncval = nc_close(ncid)))
		ERR(retncval);
	//delte 3D array	
	std::cout << "Sucess creating and storing dimension vars in: " << FileName << std::endl;
	//fflush(stdout); 
	return 0;
}
// 10.6.17 for reading point site variable arrays from nc
__device__ __host__
int read_Point_SiteVars_NC(const char* FILE_NAME, const char* varName, float* &pvar_in) //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 
											//array of dimensions
											//int pdimids; //NC_MAX_DIMS]; 1D (point) file only being read here; expected to get error message otherwise
	//const char* varNames[7] = { "slope","aspect","cc","hcan","lai","lat","lon" };
	//Open the file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, varName, &pvarid)))
		ERR(retncval);
	//read variable (input data)	
	if ((retncval = nc_get_var_float(ncid, pvarid, &pvar_in[0])))
		ERR(retncval);

	//close netcdf file	
	if ((retncval = nc_close(ncid)))
		ERR(retncval);

	return 0;
}
// 10.6.17 for reading point site variable arrays from nc
__device__ __host__
int read_Point_SiteVars_NC(const char* FILE_NAME, float** &pvar_in) //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 
											//array of dimensions
											//int pdimids; //NC_MAX_DIMS]; 1D (point) file only being read here; expected to get error message otherwise
    const char* varNames[7] = { "slope","aspect","cc","hcan","lai","lat","lon" };
											//Open the file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR(retncval);
	for (int iv = 0; iv < 7; iv++)
	{
		// get variable id
		if ((retncval = nc_inq_varid(ncid, varNames[iv], &pvarid)))
			ERR(retncval);
		//read variable (input data)	
		if ((retncval = nc_get_var_float(ncid, pvarid, &pvar_in[iv][0])))
			ERR(retncval);
	}
	//close netcdf file	
	if ((retncval = nc_close(ncid)))
		ERR(retncval);

	return 0;
}
//1.18.18 for forcing at obs points
//  read wole matrix at once
__device__ __host__
int read2DNC_Contigious(const char* FILE_NAME, const char* VAR_NAME, float** &pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR2(retncval, FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);

	//read variable (input data)	
	if (retncval = nc_get_var_float(ncid, pvarid, &pvar_in[0][0]))
		ERR(retncval);

	//close netcdf file			
	if (retncval = nc_close(ncid))
		ERR(retncval);

	return 0;
}
//10.6.17 read multiple arrays (vectors of vals at points) for given time range
__device__ __host__
int readNC_vector(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int &tStart, int tEnd, float** &pvar_in, int &nrecords)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//float* pvarin_temp = NULL;
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
											//variable data type
	nc_type varType;
	size_t pdim_sizes;
	//array of dimensions
	int pdimids[2]; // 2D file (time, point-location/index) only being read here; expected to get error message otherwise
					//dimension names 
	char pdim_Names[80];
	size_t start[2], count[2];
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR2(retncval, FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
		ERR(retncval);

	//check dimension info and set start and count arrays; 
	for (int i = 0; i < 2; i++) {
		if (retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes))
			ERR(retncval);
		if (strcmp(pdim_Names, tcor_NAME) == 0) {
			start[i] = tStart;
			if (tEnd < pdim_sizes) {
				count[i] = tEnd - tStart;
				tStart += count[i];                //new start point
			}
			else {
				count[i] = pdim_sizes - tStart;        //take the lower of the tEnd/count[i] to guarantee against going out of bounds
													   //numNc++;                               //next time go to the next netcdf file
				tStart = 0;
			}
			nrecords = count[i];
		}
		else {
			start[i] = 0;
			count[i] = pdim_sizes;
			//yxDim *= count[i];
		}
	}
	/*if (pvar_in != NULL)
	delete3DArrayblock_Contiguous(pvar_in);
	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);*/
	//read var data
	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0][0]))
		ERR(retncval);

	//close netcdf file			
	if (retncval = nc_close(ncid))
		ERR(retncval);
	return 0;
}
//2.19.18 read array at an index (time step)
__device__ __host__
int readNC_Array_atIndex(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int tIndex, float* pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//float* pvarin_temp = NULL;
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
											//variable data type
	nc_type varType;
	size_t pdim_sizes;
	//array of dimensions
	int pdimids[2]; // 2D file (time, point-location/index) only being read here; expected to get error message otherwise
					//dimension names 
	char pdim_Names[80];
	size_t start[2], count[2];
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR2(retncval, FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
		ERR(retncval);

	//check dimension info and set start and count arrays; 
	for (int i = 0; i < 2; i++) {
		if (retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes))
			ERR(retncval);
		if (strcmp(pdim_Names, tcor_NAME) == 0) {
			start[i] = tIndex;
			count[i] = 1;	
		}
		else {
			start[i] = 0;
			count[i] = pdim_sizes;
			//yxDim *= count[i];
		}
	}
	/*if (pvar_in != NULL)
	delete3DArrayblock_Contiguous(pvar_in);
	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);*/
	//read var data
	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0]))
		ERR(retncval);

	//close netcdf file			
	if (retncval = nc_close(ncid))
		ERR(retncval);
	return 0;
}
//for single point //2.19.18 read array at an index (time step)
__device__ __host__
int readNC_Array_atIndices(const char* FILE_NAME, const char* VAR_NAME, const char* tcor_NAME, int tIndex, int pIndex, float* &pvar_in)  //, MPI::Intracomm inpComm, MPI::Info inpInfo)
{
	//float* pvarin_temp = NULL;
	//ids for variable, axes,...
	int retncval = 0, ncid = 0, pvarid = 0; // pxid = 0, pyid = 0, ndims = 0; 	
											//variable data type
	nc_type varType;
	size_t pdim_sizes;
	//array of dimensions
	int pdimids[2]; // 2D file (time, point-location/index) only being read here; expected to get error message otherwise
					//dimension names 
	char pdim_Names[80];
	size_t start[2], count[2];
	//Open the netcdf file.  
	if ((retncval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
		ERR2(retncval, FILE_NAME);
	// get variable id
	if ((retncval = nc_inq_varid(ncid, VAR_NAME, &pvarid)))
		ERR(retncval);
	//var information, checking the dimension array
	if ((retncval = nc_inq_var(ncid, pvarid, NULL, &varType, NULL, pdimids, NULL)))
		ERR(retncval);

	//check dimension info and set start and count arrays; 
	for (int i = 0; i < 2; i++) {
		if (retncval = nc_inq_dim(ncid, pdimids[i], pdim_Names, &pdim_sizes))
			ERR(retncval);
		if (strcmp(pdim_Names, tcor_NAME) == 0) {
			start[i] = tIndex;
			count[i] = 1;
		}
		else {
			start[i] = pIndex;
			count[i] = 1;    // pdim_sizes;
			//yxDim *= count[i];
		}
	}
	/*if (pvar_in != NULL)
	delete3DArrayblock_Contiguous(pvar_in);
	pvar_in = create3DArrayblock_Contiguous(count[0], count[1], count[2]);*/
	//read var data
	if (retncval = nc_get_vara_float(ncid, pvarid, start, count, &pvar_in[0]))
		ERR(retncval);

	//close netcdf file			
	if (retncval = nc_close(ncid))
		ERR(retncval);
	return 0;
}


//
#include "uebDAsrPFuebpgdecls.h"
//#include "uebDAdafunctions.h"
//#include <time.h>
//#include <queue>
#pragma warning(disable : 4996)

void cuda_checkERR(cudaError_t err)
{
	if (err != cudaSuccess) {
		std::cout << "Error: " << cudaGetErrorString(err) << std::endl;
		exit(EXIT_FAILURE);
	}
}

void checkDeviceMemory()
{
	//mem check on device
	size_t freeM, totalM;
	float freeMB, totalMB, allocMB;
	cudaMemGetInfo((size_t*)&freeM, (size_t*)&totalM);
	freeMB = (size_t)freeM / (1024 * 1024);
	totalMB = (size_t)totalM / (1024 * 1024);
	allocMB = totalMB - freeMB;
	printf(" %f  MB of   %f   MB total available device memory allocated. Remaining memory =   %f MB\n", allocMB, totalMB, freeMB);
}

void estimateThroughput(size_t dataSize, clock_t beginTime, clock_t endTime)
{
	double bandWidth = dataSize * 2.0;
	double GFLOPs = (double)(bandWidth * CLOCKS_PER_SEC) / (double)(endTime - beginTime);
	printf(" Estimated throughput =  %lf GFLOPs\n", GFLOPs);
}

__global__
void callUEBrunDeviceEnsemble_HistForcing(uebCellDA *uebCellArray, int nCells, int nEns)
{
	int indx = blockIdx.x*blockDim.x + threadIdx.x;
	if (indx < nCells)
	{
		uebCellArray[indx].runUEBEnsembles_HistForcing(nEns);
	}
}
__global__
void callUEBrunDevice(uebCellDA *uebCellArray, float *uebForcArray, int nCells)
{
	int indx = blockIdx.x*blockDim.x + threadIdx.x;
	if (indx < nCells)
	{
		uebCellArray[indx].setForcingAndRadiationParamterization(uebForcArray[indx], uebForcArray[nCells + indx], uebForcArray[2 * nCells + indx],
			uebForcArray[3 * nCells + indx], uebForcArray[4 * nCells + indx], uebForcArray[5 * nCells + indx]);
		uebCellArray[indx].runUEB();
	}
}
__global__
void callUEBrunDeviceEnsemble(uebCellDA *uebCellArray, float *uebForcArray, int nCells, int nEns)
{
	int indx = blockIdx.x*blockDim.x + threadIdx.x;
	if (indx < nCells)
	{
		//cudaMalloc(&dev_loOutArr, 700000 * sizeof(float));		
		uebCellArray[indx].setForcingAndRadiationParamterization(uebForcArray[indx], uebForcArray[nCells + indx], uebForcArray[2 * nCells + indx],
			uebForcArray[3 * nCells + indx], uebForcArray[4 * nCells + indx], uebForcArray[5 * nCells + indx]);
		uebCellArray[indx].runUEBEnsembles(nEns);
	}
}
__global__
void callUEBrunDeviceEnsemble(uebCellDA *uebCellArray, float *uebForcArray, int nCells, int nEns, int nOpts,
	float *dev_multivarNormalDistSamplesForc, 
	float *dev_multivarNormalDistSamplesForcTVRH0, float *dev_multivarNormalDistSamplesForcTVRH1, float *dev_multivarNormalDistSamplesForcTVRH2)
{
	int indx = blockIdx.x*blockDim.x + threadIdx.x;
	if (indx < nCells)
	{
		//cudaMalloc(&dev_loOutArr, 700000 * sizeof(float));		
		uebCellArray[indx].copyEnsForcingMultiplier(indx, nEns, nCells + nOpts, dev_multivarNormalDistSamplesForc,
			dev_multivarNormalDistSamplesForcTVRH0, dev_multivarNormalDistSamplesForcTVRH1, dev_multivarNormalDistSamplesForcTVRH2);
		//printf(" %d forcing multi-copied ", indx);
		uebCellArray[indx].setForcingAndRadiationParamterization(uebForcArray[indx], uebForcArray[nCells + indx], uebForcArray[2 * nCells + indx],
			uebForcArray[3 * nCells + indx], uebForcArray[4 * nCells + indx], uebForcArray[5 * nCells + indx]);
		uebCellArray[indx].runUEBEnsembles(nEns);
	}
}
int callUEBrunHostEnsemble(int threadsPerBlock, int npix, int nOpts, int nEns,
	uebCellDA *uebCellDAObjArr, uebCellDA *daUebCellObjArr, float *uebForcArray, float* uebForcArrayDaPoint,
	Eigen::Matrix<float, Dynamic, Dynamic, RowMajor>  ensForcingMultiplier, std::vector<Eigen::Matrix<float, Dynamic, Dynamic, RowMajor> > ensForcingMultiplierVRH)  //,uebEnKFDA objEnKFArr      
{
	int curDv, curStr;
	//memory on device 
	//checkDeviceMemory()
	cudaError_t err = cudaSuccess;
	err = cudaGetDevice(&curDv);
	cuda_checkERR(err);
	//std::cout << " current device  = " << curDv << std::endl;
	cudaStream_t oStream;
	//cudaSetDevice(cudIndx);
	err = cudaStreamCreate(&oStream);
	cuda_checkERR(err);
	//std::cout << " current stream  = " << oStream << std::endl;
	//memory on device    
	//checkDeviceMemory();

	//float uebSizeTotMB = (float)npix*sizeof(uebCellDA) / (1024.0 * 1024.0);
	//std::cout << " Total UEB cells size in MB: " << uebSizeTotMB << std::endl;
	uebCellDA *dev_uebCellArr = NULL;
	err = cudaMalloc(&dev_uebCellArr, npix*sizeof(uebCellDA));
	cuda_checkERR(err);
	//std::cout << " device memory alloc" << std::endl;
	err = cudaMemcpyAsync(dev_uebCellArr, uebCellDAObjArr, npix*sizeof(uebCellDA), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device ueb obj copy " << std::endl;
	float *dev_uebForcArray = NULL;
	err = cudaMalloc(&dev_uebForcArray, (6 * npix)*sizeof(float));
	cuda_checkERR(err);
	//std::cout << " here 1 " << std::endl;
	err = cudaMemcpyAsync(dev_uebForcArray, uebForcArray, (6 * npix)*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " here 2 " << std::endl;
	//
	float *dev_multivarNormalDistSamplesForc = NULL;
	err = cudaMalloc(&dev_multivarNormalDistSamplesForc, (3 * (npix + nOpts) * nEns)*sizeof(float));
	cuda_checkERR(err);
	//double *host_multivarNormalDistSamplesForc = new double[3 * (npix + nOpts) * nEns];
	//host_multivarNormalDistSamplesForc = ensForcingMultiplier.data();
	err = cudaMemcpyAsync(dev_multivarNormalDistSamplesForc, ensForcingMultiplier.data(), (3 * (npix + nOpts) * nEns)*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " here 3 " << std::endl;
	//
	float *dev_multivarNormalDistSamplesForcTVRH0 = NULL;
	err = cudaMalloc(&dev_multivarNormalDistSamplesForcTVRH0, ((npix + nOpts) * nEns)*sizeof(float));
	cuda_checkERR(err);
	//double *host_ensForcingMultiplierVRH0 = new double[(npix + nOpts) * nEns];
	//host_ensForcingMultiplierVRH0 = ensForcingMultiplierVRH[0].data();
	err = cudaMemcpyAsync(dev_multivarNormalDistSamplesForcTVRH0, ensForcingMultiplierVRH[0].data(), ((npix + nOpts) * nEns)*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " here 4 " << std::endl;
	//
	float *dev_multivarNormalDistSamplesForcTVRH1 = NULL;
	err = cudaMalloc(&dev_multivarNormalDistSamplesForcTVRH1, ((npix + nOpts) * nEns)*sizeof(float));
	cuda_checkERR(err);
	//double *host_ensForcingMultiplierVRH1 = new double[(npix + nOpts) * nEns];
	//host_ensForcingMultiplierVRH1 = ensForcingMultiplierVRH[1].data();
	err = cudaMemcpyAsync(dev_multivarNormalDistSamplesForcTVRH1, ensForcingMultiplierVRH[1].data(), ((npix + nOpts) * nEns)*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " here 5 " << std::endl;
	//
	float *dev_multivarNormalDistSamplesForcTVRH2 = NULL;
	err = cudaMalloc(&dev_multivarNormalDistSamplesForcTVRH2, ((npix + nOpts) * nEns)*sizeof(float));
	cuda_checkERR(err);
	//double *host_ensForcingMultiplierVRH2 = new double[(npix + nOpts) * nEns];
	//host_ensForcingMultiplierVRH2 = ensForcingMultiplierVRH[2].data();
	err = cudaMemcpyAsync(dev_multivarNormalDistSamplesForcTVRH0, ensForcingMultiplierVRH[2].data(), ((npix + nOpts) * nEns)*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device data copy " << std::endl;
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//memory on device    checkDeviceMemory()	
	// Launch Kernel
	int blocksPerGrid = (npix + threadsPerBlock - 1) / threadsPerBlock;
	//call device run function	
	callUEBrunDeviceEnsemble << < blocksPerGrid, threadsPerBlock, 0, oStream >> > (dev_uebCellArr, dev_uebForcArray, npix, nEns, nOpts,
		dev_multivarNormalDistSamplesForc,
		dev_multivarNormalDistSamplesForcTVRH0, dev_multivarNormalDistSamplesForcTVRH1, dev_multivarNormalDistSamplesForcTVRH2);
	//synchronization
	err = cudaStreamSynchronize(oStream);   /// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//computeRun_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " finished device compute tasks for model grid" << std::endl;
	//copy data back
	err = cudaMemcpyAsync(uebCellDAObjArr, dev_uebCellArr, npix*sizeof(uebCellDA), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout << " finished device compute tasks for model grid 2" << std::endl;
	//this occurs while data copying back to host
	for (int irank = 0; irank < nOpts; irank++)
	{
		//daUebCellObjArr[irank].copyEnsForcingMultiplier(irank + npix, nEns, npix + nOpts, host_multivarNormalDistSamplesForc,
		//	host_ensForcingMultiplierVRH0, host_ensForcingMultiplierVRH1, host_ensForcingMultiplierVRH2);	
		daUebCellObjArr[irank].copyEnsForcingMultiplier(irank + npix, nEns, npix + nOpts, ensForcingMultiplier, ensForcingMultiplierVRH);
		//10.7.17 every proc. has copy of the daCellArr --- no need for broadcast
		//printf(" %d forcing multi-copied ", irank);
		daUebCellObjArr[irank].setForcingAndRadiationParamterization(uebForcArrayDaPoint[irank], uebForcArrayDaPoint[nOpts + irank], uebForcArrayDaPoint[2 * nOpts + irank],
			uebForcArrayDaPoint[3 * nOpts + irank], uebForcArrayDaPoint[4 * nOpts + irank], uebForcArrayDaPoint[5 * nOpts + irank]);
		//
		daUebCellObjArr[irank].runUEBEnsembles(nEns);

		//11.23.18 each obs points da without its observation (Leave one approach)
		daUebCellObjArr[irank + nOpts].copyEnsForcingMultiplier(irank + npix, nEns, npix + nOpts, ensForcingMultiplier, ensForcingMultiplierVRH);

		daUebCellObjArr[irank + nOpts].setForcingAndRadiationParamterization(uebForcArrayDaPoint[irank], uebForcArrayDaPoint[nOpts + irank], uebForcArrayDaPoint[2 * nOpts + irank],
			uebForcArrayDaPoint[3 * nOpts + irank], uebForcArrayDaPoint[4 * nOpts + irank], uebForcArrayDaPoint[5 * nOpts + irank]);
		//
		daUebCellObjArr[irank + nOpts].runUEBEnsembles(nEns);
	}
	//printf(" finished host compute tasks \n");
	err = cudaStreamSynchronize(oStream);// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " data copied to host" << std::endl;	//cout << err << endl;	
	/*//end of gpu call

	//EnKF
	//std::cout << "here 6 " << std::endl;
	*/
	//free device memory
	err = cudaFree(dev_uebCellArr); // != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout << " ueb obj device memory freed" << std::endl;
	err = cudaFree(dev_uebForcArray); // != cudaSuccess)	
	cuda_checkERR(err);
	//
	err = cudaFree(dev_multivarNormalDistSamplesForc); // != cudaSuccess)	
	cuda_checkERR(err);
	err = cudaFree(dev_multivarNormalDistSamplesForcTVRH0); // != cudaSuccess)	
	cuda_checkERR(err);
	err = cudaFree(dev_multivarNormalDistSamplesForcTVRH1); // != cudaSuccess)	
	cuda_checkERR(err);
	err = cudaFree(dev_multivarNormalDistSamplesForcTVRH2); // != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout <<" ueb forc device memory freed" << std::endl;

	//memory on device    
	//checkDeviceMemory();

	//clear stream 
	err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);
	//std::cout << " freed device memory" << std::endl;

	return 0;
}//
int callUEBrunHostEnsemble_HistForcing(int threadsPerBlock, int npix, int nEns, uebCellDA *uebCellDAObjArr, int nOpts, uebCellDA *daUebCellObjArr)		//, float* uebForcArrayDaPoint)  //,uebEnKFDA objEnKFArr      
{
	int curDv, curStr;
	//memory on device 
	//checkDeviceMemory()
	cudaError_t err = cudaSuccess;
	err = cudaGetDevice(&curDv);
	cuda_checkERR(err);
	//std::cout << " current device  = " << curDv << std::endl;
	cudaStream_t oStream;
	//cudaSetDevice(cudIndx);
	err = cudaStreamCreate(&oStream);
	cuda_checkERR(err);
	//std::cout << " current stream  = " << oStream << std::endl;
	//memory on device    
	//checkDeviceMemory();

	//float uebSizeTotMB = (float)npix*sizeof(uebCellDA) / (1024.0 * 1024.0);
	//std::cout << " Total UEB cells size in MB: " << uebSizeTotMB << std::endl;
	uebCellDA *dev_uebCellArr = NULL;
	err = cudaMalloc(&dev_uebCellArr, npix*sizeof(uebCellDA));
	cuda_checkERR(err);
	//std::cout << " device memory alloc" << std::endl;
	err = cudaMemcpyAsync(dev_uebCellArr, uebCellDAObjArr, npix*sizeof(uebCellDA), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device data copy " << std::endl;
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//memory on device    checkDeviceMemory()	
	// Launch Kernel
	int blocksPerGrid = (npix + threadsPerBlock - 1) / threadsPerBlock;
	//call device run function	
	callUEBrunDeviceEnsemble_HistForcing << < blocksPerGrid, threadsPerBlock, 0, oStream >> > (dev_uebCellArr, npix, nEns);
	//synchronization
	err = cudaStreamSynchronize(oStream);   // cudaDeviceSynchronize();
	cuda_checkERR(err);
	//computeRun_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " finished device compute tasks for model grid" << std::endl;
	//copy data back
	err = cudaMemcpyAsync(uebCellDAObjArr, dev_uebCellArr, npix*sizeof(uebCellDA), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
	cuda_checkERR(err);
	//this occurs while data copying back to host
	for (int irank = 0; irank < nOpts; irank++)
	{
		//9.7.18: use forc from nearest basin grid  every proc. has copy of the daCellArr --- no need for broadcast
		//daUebCellObjArr[irank].setForcingAndRadiationParamterization(uebForcArrayDaPoint[irank], uebForcArrayDaPoint[nOpts + irank], uebForcArrayDaPoint[2 * nOpts + irank],
		//	uebForcArrayDaPoint[3 * nOpts + irank], uebForcArrayDaPoint[4 * nOpts + irank], uebForcArrayDaPoint[5 * nOpts + irank]);
		//
		daUebCellObjArr[irank].runUEBEnsembles_HistForcing(nEns); //runUEBEnsembles(nEns);
		//11.23.18 each obs points da without its observation (Leave one approach)
		daUebCellObjArr[irank + nOpts].runUEBEnsembles_HistForcing(nEns);
	}
	err = cudaStreamSynchronize(oStream);// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " data copied to host" << std::endl;	//cout << err << endl;	
	/*//end of gpu call

	//EnKF
	//std::cout << "here 6 " << std::endl;
	*/
	//free device memory
	err = cudaFree(dev_uebCellArr); // != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout << " ueb obj device memory freed" << std::endl;

	//memory on device    
	//checkDeviceMemory();

	//clear stream 
	err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);

	//std::cout << " freed device memory" << std::endl;
	return 0;

}//
int callUEBrunHostEnsemble_PertForcing(int threadsPerBlock, int npix, int nEns, uebCellDA *uebCellDAObjArr, float *uebForcArray)  //,uebEnKFDA objEnKFArr      
{
	int curDv, curStr;
	//memory on device 
	//checkDeviceMemory()
	cudaError_t err = cudaSuccess;
	err = cudaGetDevice(&curDv);
	cuda_checkERR(err);
	//std::cout << " current device  = " << curDv << std::endl;
	cudaStream_t oStream;
	//cudaSetDevice(cudIndx);
	err = cudaStreamCreate(&oStream);
	cuda_checkERR(err);
	//std::cout << " current stream  = " << oStream << std::endl;
	//memory on device    
	//checkDeviceMemory();

	//float uebSizeTotMB = (float)npix*sizeof(uebCellDA) / (1024.0 * 1024.0);
	//std::cout << " Total UEB cells size in MB: " << uebSizeTotMB << std::endl;
	uebCellDA *dev_uebCellArr = NULL;
	err = cudaMalloc(&dev_uebCellArr, npix*sizeof(uebCellDA));
	cuda_checkERR(err);
	//std::cout << " device memory alloc" << std::endl;
	err = cudaMemcpyAsync(dev_uebCellArr, uebCellDAObjArr, npix*sizeof(uebCellDA), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device ueb obj copy " << std::endl;
	float *dev_uebForcArray = NULL;
	err = cudaMalloc(&dev_uebForcArray, (6 * npix)*sizeof(float));
	cuda_checkERR(err);
	err = cudaMemcpyAsync(dev_uebForcArray, uebForcArray, (6 * npix)*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device data copy " << std::endl;
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//memory on device    checkDeviceMemory()	
	// Launch Kernel
	int blocksPerGrid = (npix + threadsPerBlock - 1) / threadsPerBlock;
	//call device run function	
	callUEBrunDeviceEnsemble << < blocksPerGrid, threadsPerBlock, 0, oStream >> > (dev_uebCellArr, dev_uebForcArray, npix, nEns);
	//synchronization
	err = cudaStreamSynchronize(oStream);   /// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//computeRun_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " finished device compute tasks for model grid" << std::endl;
	//copy data back
	err = cudaMemcpyAsync(uebCellDAObjArr, dev_uebCellArr, npix*sizeof(uebCellDA), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
	cuda_checkERR(err);

	err = cudaStreamSynchronize(oStream);// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " data copied to host" << std::endl;	//cout << err << endl;	
	/*//end of gpu call

	//EnKF
	//std::cout << "here 6 " << std::endl;
	*/
	//free device memory
	err = cudaFree(dev_uebCellArr); // != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout << " ueb obj device memory freed" << std::endl;
	err = cudaFree(dev_uebForcArray); // != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout <<" ueb forc device memory freed" << std::endl;

	//memory on device    
	//checkDeviceMemory();

	//clear stream 
	err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);

	//std::cout << " freed device memory" << std::endl;

	return 0;
}//
int callUEBrunHost(int threadsPerBlock, int npix, uebCellDA *uebCellDAObjArr, float *uebForcArray)  //,uebEnKFDA objEnKFArr      
{
	int curDv, curStr;
	//memory on device 
	//checkDeviceMemory()
	cudaError_t err = cudaSuccess;
	err = cudaGetDevice(&curDv);
	cuda_checkERR(err);
	//std::cout << " current device  = " << curDv << std::endl;
	cudaStream_t oStream;
	//cudaSetDevice(cudIndx);
	err = cudaStreamCreate(&oStream);
	cuda_checkERR(err);
	//std::cout << " current stream  = " << oStream << std::endl;
	//memory on device    
	//checkDeviceMemory();

	//float uebSizeTotMB = (float)npix*sizeof(uebCellDA) / (1024.0 * 1024.0);
	//std::cout << " Total UEB cells size in MB: " << uebSizeTotMB << std::endl;

	uebCellDA *dev_uebCellArr = NULL;
	err = cudaMalloc(&dev_uebCellArr, npix*sizeof(uebCellDA));
	cuda_checkERR(err);
	//std::cout << " device memory alloc" << std::endl;
	err = cudaMemcpyAsync(dev_uebCellArr, uebCellDAObjArr, npix*sizeof(uebCellDA), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device ueb obj copy " << std::endl;
	float *dev_uebForcArray = NULL;
	err = cudaMalloc(&dev_uebForcArray, (6 * npix)*sizeof(float));
	cuda_checkERR(err);
	err = cudaMemcpyAsync(dev_uebForcArray, uebForcArray, (6 * npix)*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();

	//std::cout << " device data copied " << std::endl;
	//memory on device    
	//checkDeviceMemory();
	// Launch Kernel
	int blocksPerGrid = (npix + threadsPerBlock - 1) / threadsPerBlock;
	//call device run function	
	callUEBrunDevice << < blocksPerGrid, threadsPerBlock, 0, oStream >> > (dev_uebCellArr, dev_uebForcArray, npix);
	//synchronization
	err = cudaStreamSynchronize(oStream);   /// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//computeRun_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();

	//std::cout << " finished device compute tasks for model grid" << std::endl;
	//copy data back
	err = cudaMemcpyAsync(uebCellDAObjArr, dev_uebCellArr, npix*sizeof(uebCellDA), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();

	//std::cout << " data copied to host" << std::endl;	//cout << err << endl;	
	/*end of gpu call
	//EnKF
	//std::cout << "here 6 " << std::endl;
	*/
	//free device memory
	err = cudaFree(dev_uebCellArr); // != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout << " ueb obj device memory freed" << std::endl;
	err = cudaFree(dev_uebForcArray); // != cudaSuccess)	
	cuda_checkERR(err);
	//std::cout << " ueb forc device memory freed" << std::endl;
	//clear stream 
	err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);

	//std::cout << " freed device memory " << std::endl;

	return 0;
}//
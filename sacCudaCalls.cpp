//
#include "uebDAsrPFuebpgdecls.h"
//#include "uebDAdafunctions.h"
//#include <time.h>
//#include <queue>
#pragma warning(disable : 4996)

void cuda_checkERR(cudaError_t err);
void checkDeviceMemory();
void estimateThroughput(size_t dataSize, clock_t beginTime, clock_t endTime);

__host__ __device__
void SACGridDA::runSACGrid_Ens(int grdIndx, int nEns, float *pxvEns)
{
	float edmndOrg = edmnd;
	for (int ie = 0; ie < nEns; ie++)
	{
		edmnd = edmndOrg;
		fland1Device(pxvEns[grdIndx * nEns + ie], dtday, aesc, sacstEns[ie], sacpar,
			edmnd, surfEns[ie], grndEns[ie], tet,
			sif, bfp, bfs, ssur, sdro, sperc, hrapx, hrapy, error);
	}
}

__global__
void do_sac_EnsDevice_Ens(int npix, int nEns, SACGridDA *sacGridArray, float *pxvEns)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < npix)
	{
		sacGridArray[i].runSACGrid_Ens(i, nEns, pxvEns);

	}  /* End of grid loop */

}   /* End of program */


int do_sac_EnsHost_Ens(int threadsPerBlock, int npix, int nEns, SACGridDA *sacGridArray, float *pxvEns)
{
	//std::cout << " here in host call_8 " << std::endl;
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

	SACGridDA *dev_sacGridArray = NULL;
	err = cudaMalloc(&dev_sacGridArray, npix*sizeof(SACGridDA));
	cuda_checkERR(err);
	//std::cout << " device memory alloc ";
	err = cudaMemcpyAsync(dev_sacGridArray, sacGridArray, npix*sizeof(SACGridDA), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//
	float *dev_pxvEns = NULL;
	err = cudaMalloc(&dev_pxvEns, nEns*npix*sizeof(float));
	cuda_checkERR(err);
	err = cudaMemcpyAsync(dev_pxvEns, pxvEns, nEns*npix*sizeof(float), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device data copy ";
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();

	//memory on device    
	//checkDeviceMemory();

	// Launch Kernel
	int blocksPerGrid = (npix + threadsPerBlock - 1) / threadsPerBlock;
	do_sac_EnsDevice_Ens << < blocksPerGrid, threadsPerBlock, 0 >> > (npix, nEns, dev_sacGridArray, dev_pxvEns);
	//synchronization
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//computeRun_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " finished device compute tasks " << std::endl;
	//copy data back
	err = cudaMemcpyAsync(sacGridArray, dev_sacGridArray, npix*sizeof(SACGridDA), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " data copied to host ";	//cout << err << endl;	
	//checkDeviceMemory();

	//free device memory
	err = cudaFree(dev_sacGridArray);                      // npsac_def*npix*sizeof(float));
	cuda_checkERR(err);
	err = cudaFree(dev_pxvEns);                      // npsac_def*npix*sizeof(float));
	cuda_checkERR(err);
	//std::cout <<" Device memory freed ";
	//memory on device    
	//checkDeviceMemory();

	//clear stream 
	err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);

	return 0;
}
//
__global__
void do_sac_EnsDevice(int npix, SACGridDA *sacGridArray)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i < npix)
	{
		sacGridArray[i].runSACGrid();
	}  /* End of grid loop */

}   /* End of program */

int do_sac_EnsHost(int threadsPerBlock, int npix, SACGridDA *sacGridArray)
{
	//std::cout << " here in host call_8 " << std::endl;
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

	SACGridDA *dev_sacGridArray = NULL;
	err = cudaMalloc(&dev_sacGridArray, npix*sizeof(SACGridDA));
	cuda_checkERR(err);
	//std::cout << " device memory alloc ";
	err = cudaMemcpyAsync(dev_sacGridArray, sacGridArray, npix*sizeof(SACGridDA), cudaMemcpyHostToDevice, oStream); // != cudaSuccess)		
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//std::cout << " device data copy ";
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();

	//memory on device    
	//checkDeviceMemory();

	// Launch Kernel
	int blocksPerGrid = (npix + threadsPerBlock - 1) / threadsPerBlock;
	do_sac_EnsDevice << < blocksPerGrid, threadsPerBlock, 0 >> > (npix, dev_sacGridArray);
	//synchronization
	err = cudaStreamSynchronize(oStream);
	cuda_checkERR(err);
	//computeRun_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " finished device compute tasks " << std::endl;
	//copy data back
	err = cudaMemcpyAsync(sacGridArray, dev_sacGridArray, npix*sizeof(SACGridDA), cudaMemcpyDeviceToHost, oStream);// != cudaSuccess)	
	cuda_checkERR(err);
	err = cudaStreamSynchronize(oStream);// cudaDeviceSynchronize();
	cuda_checkERR(err);
	//dataCopy_Time += (MPI::Wtime() - intermStart_Time);
	//intermStart_Time = MPI::Wtime();
	//std::cout << " data copied to host ";	//cout << err << endl;	
	//checkDeviceMemory();

	//free device memory
	err = cudaFree(dev_sacGridArray);                      // npsac_def*npix*sizeof(float));
	cuda_checkERR(err);
	//std::cout <<" Device memory freed ";
	//memory on device    
	//checkDeviceMemory();

	//clear stream 
	err = cudaStreamDestroy(oStream);
	cuda_checkERR(err);

	return 0;
}

void do_sac_EnsHost_PC(int npix, SACGridDA *sacGridArray)
{
	for (int i = 0; i<npix; i++)
		//if (i < npix)
	{
		sacGridArray[i].runSACGrid();
	}  /* End of grid loop */

}   /* End of program */
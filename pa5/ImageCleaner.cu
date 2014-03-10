#include "ImageCleaner.h"

#ifndef SIZEX
#define SIZEX    1024
#endif
#ifndef SIZEY
#define SIZEY    1024
#endif
#define BLOCK_SIZE 16
#define PI 3.14159265

typedef struct {
	int size;
	int stride;
	float *elements;
} Matrix;
//----------------------------------------------------------------
// TODO:  CREATE NEW KERNELS HERE.  YOU CAN PLACE YOUR CALLS TO
//        THEM IN THE INDICATED SECTION INSIDE THE 'filterImage'
//        FUNCTION.
//
// BEGIN ADD KERNEL DEFINITIONS
//----------------------------------------------------------------

//row and col is within the sub-matrix
__device__ float getElement(const Matrix A, int row, int col) {
	return A.elements[row * A.stride + col];
}

__device__ void setElement(Matrix A, int row, int col, float value) {
	A.elements[row * A.stride + col] = value;
}

__device__ Matrix getSubMatrix(Matrix A, int blockRow, int blockCol) {
	Matrix Asub;
	Asub.size = BLOCK_SIZE;
	Asub.stride = A.stride;
	Asub.elements = A.elements + blockRow * BLOCK_SIZE * A.stride + blockCol * BLOCK_SIZE;
	return Asub;
}

__global__ void FT_Kernel(const Matrix realMatrix, const Matrix imagMatrix, const Matrix sinMatrix,
		const Matrix cosMatrix, Matrix resultRealMatrix, Matrix resultImagMatrix) {
	int blockRow = blockIdx.y;
	int blockCol = blockIdx.x;
	// the sub-matrix block this thread block is responsible for
	Matrix subResultReal = getSubMatrix(resultRealMatrix, blockCol, blockRow);
	Matrix subResultImag = getSubMatrix(resultImagMatrix, blockCol, blockRow);
	//the sub-matrix entry this thread is responsible for
	int row = threadIdx.y;
	int col = threadIdx.x;

	float Creal = 0;
	float Cimag = 0;

	for (int m = 0; m < (realMatrix.size / BLOCK_SIZE); ++m) {
		Matrix realSub = getSubMatrix(realMatrix, blockRow, m);
		Matrix imagSub = getSubMatrix(imagMatrix, blockRow, m);
		Matrix sinSub = getSubMatrix(sinMatrix, m, blockCol);
		Matrix cosSub = getSubMatrix(cosMatrix, m, blockCol);

		__shared__
		float realShare[BLOCK_SIZE][BLOCK_SIZE];
		__shared__
		float imagShare[BLOCK_SIZE][BLOCK_SIZE];
		__shared__
		float sinShare[BLOCK_SIZE][BLOCK_SIZE];
		__shared__
		float cosShare[BLOCK_SIZE][BLOCK_SIZE];

		realShare[row][col] = getElement(realSub, row, col);
		imagShare[row][col] = getElement(imagSub, row, col);
		sinShare[row][col] = getElement(sinSub, row, col);
		cosShare[row][col] = getElement(cosSub, row, col);

		__syncthreads();

		for (int i = 0; i < BLOCK_SIZE; ++i) {
			Creal += realShare[row][i] * cosShare[i][col] + imagShare[row][i] * sinShare[i][col];
			Cimag += realShare[row][i] * sinShare[i][col]*(-1) + imagShare[row][i] * cosShare[i][col];
		}
		__syncthreads();
	}
	//Write sub-result to device memory
	//each thread writes one element
	//must write the result to a buffer instead of write back, since other blocks might have not finished.
	setElement(subResultReal, col, row, Creal);
	setElement(subResultImag, col, row, Cimag);

}

__global__ void iFT_Kernel(const Matrix realMatrix, const Matrix imagMatrix, const Matrix sinMatrix,
		const Matrix cosMatrix, Matrix resultRealMatrix, Matrix resultImagMatrix) {
	//blockIdx and threadIdx are there as keyword
	int blockRow = blockIdx.y;
	int blockCol = blockIdx.x;
	// the sub-matrix block this thread block is responsible for
	Matrix subResultReal = getSubMatrix(resultRealMatrix, blockCol, blockRow);
	Matrix subResultImag = getSubMatrix(resultImagMatrix, blockCol, blockRow);
	//the sub-matrix entry this thread is responsible for
	int row = threadIdx.y;
	int col = threadIdx.x;

	float Creal = 0;
	float Cimag = 0;

	for (int m = 0; m < (realMatrix.size / BLOCK_SIZE); ++m) {
		Matrix realSub = getSubMatrix(realMatrix, blockRow, m);
		Matrix imagSub = getSubMatrix(imagMatrix, blockRow, m);
		Matrix sinSub = getSubMatrix(sinMatrix, m, blockCol);
		Matrix cosSub = getSubMatrix(cosMatrix, m, blockCol);

		__shared__
		float realShare[BLOCK_SIZE][BLOCK_SIZE];
		__shared__
		float imagShare[BLOCK_SIZE][BLOCK_SIZE];
		__shared__
		float sinShare[BLOCK_SIZE][BLOCK_SIZE];
		__shared__
		float cosShare[BLOCK_SIZE][BLOCK_SIZE];

		realShare[row][col] = getElement(realSub, row, col);
		imagShare[row][col] = getElement(imagSub, row, col);
		sinShare[row][col] = getElement(sinSub, row, col);
		cosShare[row][col] = getElement(cosSub, row, col);

		__syncthreads();

		for (int i = 0; i < BLOCK_SIZE; ++i) {
			Creal += realShare[row][i] * cosShare[i][col] - imagShare[row][i] * sinShare[i][col];
			Cimag += imagShare[row][i] * cosShare[i][col] + realShare[row][i] * sinShare[i][col];
		}
		__syncthreads();
	}
	//Write sub-result to device memory
	//each thread writes one element
	setElement(subResultReal, col, row, Creal / SIZEX);
	setElement(subResultImag, col, row, Cimag / SIZEX);

}

__global__ void sincosKernel(Matrix sinMatrix, Matrix cosMatrix) {
	int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
	int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
	float angel = row * col * 2 * 3.1415926536 / sinMatrix.size;
	//each thread is responsible for filling in only one entry
	sinMatrix.elements[row * sinMatrix.size + col] = __sinf(angel);
	cosMatrix.elements[row * cosMatrix.size + col] = __cosf(angel);

}

__global__ void filterKernel(Matrix realMatrix, Matrix imagMatrix) {
	int row = blockIdx.y * BLOCK_SIZE + threadIdx.y;
	int col = blockIdx.x * BLOCK_SIZE + threadIdx.x;
	int low8 = realMatrix.size / 8;
	int high8 = realMatrix.size - low8;
	if (!((row < low8 && col < low8) || (row < low8 && col >= high8) || (row >= high8 && col >= high8)
			|| (row >= high8 && col < low8))) {
		realMatrix.elements[row * realMatrix.size + col] = 0;
		imagMatrix.elements[row * imagMatrix.size + col] = 0;
	}
}

//----------------------------------------------------------------
// END ADD KERNEL DEFINTIONS
//----------------------------------------------------------------

__host__ float filterImage(float *real_image, float *imag_image, int size_x, int size_y) {
// check that the sizes match up
	assert(size_x == SIZEX);
	assert(size_y == SIZEY);

	int matSize = size_x * size_y * sizeof(float);

// These variables are for timing purposes
	float transferDown = 0, transferUp = 0, execution = 0;
	cudaEvent_t start, stop;

	CUDA_ERROR_CHECK(cudaEventCreate(&start));
	CUDA_ERROR_CHECK(cudaEventCreate(&stop));

// Create a stream and initialize it
	cudaStream_t filterStream;
	CUDA_ERROR_CHECK(cudaStreamCreate(&filterStream));

// Alloc space on the device
	Matrix realMatrix, imagMatrix, tmpRealMatrix, tmpImagMatrix, sinMatrix, cosMatrix;
	realMatrix.size = imagMatrix.size = tmpRealMatrix.size = tmpImagMatrix.size = sinMatrix.size = cosMatrix.size =
			size_x;
	realMatrix.stride = imagMatrix.stride = tmpRealMatrix.stride = tmpImagMatrix.stride = sinMatrix.stride =
			cosMatrix.stride = size_x;
	//device memory could be allocated in the host, but not dereference in the host
	CUDA_ERROR_CHECK(cudaMalloc((void** )&realMatrix.elements, matSize));
	CUDA_ERROR_CHECK(cudaMalloc((void** )&imagMatrix.elements, matSize));

	CUDA_ERROR_CHECK(cudaMalloc((void** )&tmpRealMatrix.elements, matSize));
	CUDA_ERROR_CHECK(cudaMalloc((void** )&tmpImagMatrix.elements, matSize));

	CUDA_ERROR_CHECK(cudaMalloc((void** )&sinMatrix.elements, matSize));
	CUDA_ERROR_CHECK(cudaMalloc((void** )&cosMatrix.elements, matSize));

// Start timing for transfer down
	CUDA_ERROR_CHECK(cudaEventRecord(start, filterStream));

// Here is where we copy matrices down to the device
	CUDA_ERROR_CHECK(cudaMemcpy(realMatrix.elements, real_image, matSize, cudaMemcpyHostToDevice));
	CUDA_ERROR_CHECK(cudaMemcpy(imagMatrix.elements, imag_image, matSize, cudaMemcpyHostToDevice));

// Stop timing for transfer down
	CUDA_ERROR_CHECK(cudaEventRecord(stop, filterStream));
	CUDA_ERROR_CHECK(cudaEventSynchronize(stop));
	CUDA_ERROR_CHECK(cudaEventElapsedTime(&transferDown, start, stop));

// Start timing for the execution
	CUDA_ERROR_CHECK(cudaEventRecord(start, filterStream));

//----------------------------------------------------------------
// TODO: YOU SHOULD PLACE ALL YOUR KERNEL EXECUTIONS
//        HERE BETWEEN THE CALLS FOR STARTING AND
//        FINISHING TIMING FOR THE EXECUTION PHASE
// BEGIN ADD KERNEL CALLS
//----------------------------------------------------------------

// This is an example kernel call, you should feel free to create
// as many kernel calls as you feel are needed for your program
// Each of the parameters are as follows:
//    1. Number of thread blocks, can be either int or dim3 (see CUDA manual)
//    2. Number of threads per thread block, can be either int or dim3 (see CUDA manual)
//    3. Always should be '0' unless you read the CUDA manual and learn about dynamically allocating shared memory
//    4. Stream to execute kernel on, should always be 'filterStream'
//
// Also note that you pass the pointers to the device memory to the kernel call
	dim3 dimGrid(size_x / BLOCK_SIZE, size_y / BLOCK_SIZE);
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);

sincosKernel<<<dimGrid, dimBlock, 0, filterStream>>> (sinMatrix, cosMatrix);
	FT_Kernel<<<dimGrid, dimBlock, 0, filterStream>>> (realMatrix, imagMatrix, sinMatrix, cosMatrix, tmpRealMatrix, tmpImagMatrix);
	FT_Kernel<<<dimGrid, dimBlock, 0, filterStream>>> (tmpRealMatrix, tmpImagMatrix, sinMatrix, cosMatrix, realMatrix, imagMatrix);
	filterKernel<<<dimGrid, dimBlock, 0, filterStream>>> (realMatrix, imagMatrix);
	iFT_Kernel<<<dimGrid, dimBlock, 0, filterStream>>> (realMatrix, imagMatrix, sinMatrix, cosMatrix, tmpRealMatrix, tmpImagMatrix);
	iFT_Kernel<<<dimGrid, dimBlock, 0, filterStream>>> (tmpRealMatrix, tmpImagMatrix, sinMatrix, cosMatrix, realMatrix, imagMatrix);


  //---------------------------------------------------------------- 
  // END ADD KERNEL CALLS
  //----------------------------------------------------------------

  // Finish timimg for the execution 
	CUDA_ERROR_CHECK(cudaEventRecord(stop, filterStream));
	CUDA_ERROR_CHECK(cudaEventSynchronize(stop));
	CUDA_ERROR_CHECK(cudaEventElapsedTime(&execution, start, stop));

// Start timing for the transfer up
	CUDA_ERROR_CHECK(cudaEventRecord(start, filterStream));

// Here is where we copy matrices back from the device
	CUDA_ERROR_CHECK(cudaMemcpy(real_image, realMatrix.elements, matSize, cudaMemcpyDeviceToHost));
	CUDA_ERROR_CHECK(cudaMemcpy(imag_image, imagMatrix.elements, matSize, cudaMemcpyDeviceToHost));

// Finish timing for transfer up
	CUDA_ERROR_CHECK(cudaEventRecord(stop, filterStream));
	CUDA_ERROR_CHECK(cudaEventSynchronize(stop));
	CUDA_ERROR_CHECK(cudaEventElapsedTime(&transferUp, start, stop));

// Synchronize the stream
	CUDA_ERROR_CHECK(cudaStreamSynchronize(filterStream));
// Destroy the stream
	CUDA_ERROR_CHECK(cudaStreamDestroy(filterStream));
// Destroy the events
	CUDA_ERROR_CHECK(cudaEventDestroy(start));
	CUDA_ERROR_CHECK(cudaEventDestroy(stop));

// Free the memory
	CUDA_ERROR_CHECK(cudaFree(realMatrix.elements));
	CUDA_ERROR_CHECK(cudaFree(imagMatrix.elements));
	CUDA_ERROR_CHECK(cudaFree(tmpRealMatrix.elements));
	CUDA_ERROR_CHECK(cudaFree(tmpImagMatrix.elements));
	CUDA_ERROR_CHECK(cudaFree(sinMatrix.elements));
	CUDA_ERROR_CHECK(cudaFree(cosMatrix.elements));

// Dump some usage statistics
	printf("CUDA IMPLEMENTATION STATISTICS:\n");
	printf("  Host to Device Transfer Time: %f ms\n", transferDown);
	printf("  Kernel(s) Execution Time: %f ms\n", execution);
	printf("  Device to Host Transfer Time: %f ms\n", transferUp);
	float totalTime = transferDown + execution + transferUp;
	printf("  Total CUDA Execution Time: %f ms\n\n", totalTime);
// Return the total time to transfer and execute
	return totalTime;
}


#include "ImageCleaner.h"
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <omp.h>
#include <xmmintrin.h>

#define PI 3.14159265
#define CHUNK_SIZE 20

void precalculate(float *sin1, float *sin3, float *cos1, float *cos3, int size_x, int size_y) {
	float term;
	int n, y, x;
//#pragma omp parallel sections
	{
//#pragma omp section
		{
#pragma omp parallel for shared(sin1, cos1, size_y, size_x) private(y, n, term) collapse(2)
			for (y = 0; y < size_y; y++) {
				for (n = 0; n < size_y; n++) {
					term = -2 * PI * y * n / size_y;
					sin1[y * size_y + n] = sin(term);
					cos1[y * size_y + n] = cos(term);
				}
			}
		}
//#pragma omp section
		{
#pragma omp parallel for shared(sin1, cos1, size_y, size_x) private(x, n, term) collapse(2)
			for (x = 0; x < size_x; x++) {
				for (n = 0; n < size_y; n++) {
					term = -2 * PI * x * n / size_x;
					sin3[x * size_x + n] = sin(term);
					cos3[x * size_x + n] = cos(term);
				}
			}
		}
	}
}



static inline void calc_out_buffer(int size_y, int y, const float *cos1, const float
		  *sin1, const float *real_image_cp, const float *imag_image_cp, 
		  float *reduction_RealOutBuffer, float *reduction_ImagOutBuffer){

  float a_RealOutBuffer[4], a_ImagOutBuffer[4];

  __m128 v_fft_real, v_fft_image, v_real_image,
	v_image_image, v_a, v_b;
  __m128 v_real_out  = _mm_set_ps1(0.f);
  __m128 v_image_out = _mm_set_ps1(0.f);

  for (int n = 0; n < size_y; n+=128/32) {
	v_fft_real    = _mm_load_ps(cos1 + n);
	v_fft_image   = _mm_load_ps(sin1 + n);
	v_real_image  = _mm_load_ps(real_image_cp + n);
	v_image_image = _mm_load_ps(imag_image_cp + n);
	//reduction_RealOutBuffer
	v_a        = _mm_mul_ps(v_real_image , v_fft_real ); 
	v_b        = _mm_mul_ps(v_image_image, v_fft_image);
	v_real_out = _mm_add_ps(v_real_out,v_a);
	v_real_out = _mm_sub_ps(v_real_out,v_b); 
	//reduction_ImagOutBuffer
	v_a         = _mm_mul_ps(v_image_image, v_fft_real); 
	v_b         = _mm_mul_ps(v_real_image, v_fft_image);
	v_image_out = _mm_add_ps(v_image_out,v_a);
	v_image_out = _mm_add_ps(v_image_out,v_b); 


	//fft_real = cos1[y * size_y + n];
	//fft_imag = sin1[y * size_y + n];
	//reduction_RealOutBuffer += (real_image_cp[n] * fft_real) - (imag_image_cp[n] * fft_imag);
	//reduction_ImagOutBuffer += (imag_image_cp[n] * fft_real) + (real_image_cp[n] * fft_imag);
  }
  _mm_store_ps(a_RealOutBuffer, v_real_out);
  _mm_store_ps(a_ImagOutBuffer, v_image_out);

 for(int m = 1;m<4;m++){
	a_RealOutBuffer[0] += a_RealOutBuffer[m];
	a_ImagOutBuffer[0] += a_ImagOutBuffer[m];
  }

  *reduction_RealOutBuffer = a_RealOutBuffer[0];
  *reduction_ImagOutBuffer = a_ImagOutBuffer[0];
 
}





//#pragma omp parallel for shared(sin1, cos1, size_y, size_x) private(y, n, term) collapse(2)

//#pragma omp parallel for shared(sin1, cos1, size_y, size_x) private(x, n, term) collapse(2)

//111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
//111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
void cpu_fftx(float *real_image, float *imag_image, int size_x, int size_y, float *sin1, float *sin3, float *cos1,
		float *cos3) {
// Create some space for storing temporary values
	float *realOutBuffer = new float[size_x];
	float *imagOutBuffer = new float[size_x];
// Local values
	float fft_real;
	float fft_imag;
	unsigned int x, y;
	//float term;
	//float *real_image_cp = new float[size_y];
	//float *imag_image_cp = new float[size_y];

	float reduction_RealOutBuffer;
	float reduction_ImagOutBuffer;

	int chunk = 200;

	int eightY = size_y / 8;
	int eight7Y = size_y - eightY;

//#pragma omp parallel for shared(real_image, imag_image, size_x, size_y, \
  //								realOutBuffer, imagOutBuffer )		\
  //private(x, y, n, term, fft_real, fft_imag,reduction_RealOutBuffer,reduction_ImagOutBuffer )
	for (x = 0; x < size_x; x++) {
	/*	for (m = 0; m < size_y; m++) {
			real_image_cp[m] = real_image[x * size_x + m];
			imag_image_cp[m] = imag_image[x * size_x + m];
		}
*/
#pragma omp parallel for shared(real_image, imag_image, size_x, size_y, x, \
		realOutBuffer, imagOutBuffer, eightY, eight7Y) \
        private(y, reduction_RealOutBuffer, \
        reduction_ImagOutBuffer) schedule(dynamic, 4)

		for (y = 0; y < size_y; y++) {
			if (!(y >= eightY && y < eight7Y)) {
calc_out_buffer(size_y, y, cos1+y * size_y, sin1+y * size_y, real_image + x * size_x, imag_image + x * size_x, &reduction_RealOutBuffer, &reduction_ImagOutBuffer);
/*
				reduction_RealOutBuffer = 0.0f;
				reduction_ImagOutBuffer = 0.0f;

				for (n = 0; n < size_y; n++) {
					term = -2 * PI * y * n / size_y;
//				fft_real = cos(term);
//				fft_imag = sin(term);

					fft_real = cos1[y * size_y + n];
					fft_imag = sin1[y * size_y + n];
					reduction_RealOutBuffer += (real_image_cp[n] * fft_real) - (imag_image_cp[n] * fft_imag);
					reduction_ImagOutBuffer += (imag_image_cp[n] * fft_real) + (real_image_cp[n] * fft_imag);
				}
*/
				realOutBuffer[y] = reduction_RealOutBuffer;
				imagOutBuffer[y] = reduction_ImagOutBuffer;
			}
			else {
				realOutBuffer[y] = 0;
				imagOutBuffer[y] = 0;
			}
		}

// Write the buffer back to were the original values were
		for (int yy = 0; yy < size_y; yy++) {
			real_image[x * size_x + yy] = realOutBuffer[yy];
			imag_image[x * size_x + yy] = imagOutBuffer[yy];
		}
	}

// Reclaim some memory
	delete[] realOutBuffer;
	delete[] imagOutBuffer;
	//delete[] real_image_cp;
	//delete[] imag_image_cp;
}

//222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
//222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
// This is the same as the thing above, except it has a scaling factor added to it
void cpu_ifftx(float *real_image, float *imag_image, int size_x, int size_y, float*sin1, float*sin3, float *cos1,
		float *cos3) {
// Create some space for storing temporary values
	float *realOutBuffer = new float[size_x];
	float *imagOutBuffer = new float[size_x];
	float fft_real;
	float fft_imag;
	unsigned int n, m, x, y;
	float term;

	float *real_image_cp = new float[size_y];
	float *imag_image_cp = new float[size_y];

	float reduction_RealOutBuffer;
	float reduction_ImagOutBuffer;

	int eightY = size_y / 8;
	int eight7Y = size_y - eightY;
	int eightX = size_x / 8;
	int eight7X = size_x - eightX;

	for (x = 0; x < size_x; x++) {
//		if (x >= eightX && x < eight7X)
//			break;

		for (m = 0; m < size_y; m++) {
			real_image_cp[m] = real_image[x * size_x + m];
			imag_image_cp[m] = imag_image[x * size_x + m];
		}

#pragma omp parallel for shared(real_image, imag_image, size_x, size_y, x, \
		realOutBuffer, imagOutBuffer, real_image_cp, imag_image_cp) \
        private(y, n, term, fft_real, fft_imag,reduction_RealOutBuffer, \
        reduction_ImagOutBuffer) schedule(dynamic, 4)

		for (y = 0; y < size_y; y++) {
			reduction_RealOutBuffer = 0.0f;
			reduction_ImagOutBuffer = 0.0f;

			for (n = 0; n < size_y; n++) {
				if (!(n >= eightY && n < eight7Y)) {
					// Compute the frequencies for this index
					float term = 2 * PI * y * n / size_y;
//				fft_real = cos(term);
//				fft_imag = sin(term);

					fft_real = cos1[y * size_y + n];
					fft_imag = -sin1[y * size_y + n];

					reduction_RealOutBuffer += (real_image_cp[n] * fft_real) - (imag_image_cp[n] * fft_imag);
					reduction_ImagOutBuffer += (imag_image_cp[n] * fft_real) + (real_image_cp[n] * fft_imag);
				}
			}
			// Incoporate the scaling factor here
			realOutBuffer[y] = reduction_RealOutBuffer / size_y;
			imagOutBuffer[y] = reduction_ImagOutBuffer / size_y;
		}

// Write the buffer back to were the original values were
		for (int yy = 0; yy < size_y; yy++) {
			real_image[x * size_x + yy] = realOutBuffer[yy];
			imag_image[x * size_x + yy] = imagOutBuffer[yy];
		}
	}
// Reclaim some memory
	delete[] realOutBuffer;
	delete[] imagOutBuffer;
	delete[] real_image_cp;
	delete[] imag_image_cp;
}

//33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
//33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
void cpu_ffty(float *real_image, float *imag_image, int size_x, int size_y, float*sin1, float*sin3, float *cos1,
		float *cos3) {
// Allocate some space for temporary values
	float *realOutBuffer = new float[size_y];
	float *imagOutBuffer = new float[size_y];
	float fft_real;
	float fft_imag;
	unsigned int n, m, x, y;
	float term;
	float *real_image_cp = new float[size_y];
	float *imag_image_cp = new float[size_y];

	float reduction_RealOutBuffer;
	float reduction_ImagOutBuffer;
	int eightX = size_x / 8;
	int eight7X = size_x - eightX;

	for (y = 0; y < size_y; y++) {

		for (m = 0; m < size_y; m++) {
			real_image_cp[m] = real_image[m * size_x + y];
			imag_image_cp[m] = imag_image[m * size_x + y];
		}

#pragma omp parallel for shared(real_image, imag_image, size_x, size_y, y, \
		realOutBuffer, imagOutBuffer, real_image_cp, imag_image_cp) \
        private(x, reduction_RealOutBuffer, \
        reduction_ImagOutBuffer) schedule(dynamic, 4)

		for (x = 0; x < size_x; x++) {
			if (!(x >= eightX && x < eight7X)) {
calc_out_buffer(size_y, y, cos3+x * size_x, sin3+x * size_x, real_image_cp, imag_image_cp, &reduction_RealOutBuffer, &reduction_ImagOutBuffer);
/*
				reduction_RealOutBuffer = 0.0f;
				reduction_ImagOutBuffer = 0.0f;

				for (n = 0; n < size_y; n++) {
					term = -2 * PI * x * n / size_x;
//				fft_real = cos(term);
//				fft_imag = sin(term);

					fft_real = cos3[x * size_x + n];
					fft_imag = sin3[x * size_x + n];

					reduction_RealOutBuffer += (real_image_cp[n] * fft_real) - (imag_image_cp[n] * fft_imag);
					reduction_ImagOutBuffer += (imag_image_cp[n] * fft_real) + (real_image_cp[n] * fft_imag);
				}
*/
				realOutBuffer[x] = reduction_RealOutBuffer;
				imagOutBuffer[x] = reduction_ImagOutBuffer;
			}
			else {
				realOutBuffer[x] = 0;
				imagOutBuffer[x] = 0;
			}
		}
// Write the buffer back to were the original values were
		for (int xx = 0; xx < size_x; xx++) {
			real_image[xx * size_x + y] = realOutBuffer[xx];
			imag_image[xx * size_x + y] = imagOutBuffer[xx];
		}
	}
// Reclaim some memory
	delete[] realOutBuffer;
	delete[] imagOutBuffer;
	delete[] real_image_cp;
	delete[] imag_image_cp;
}

//44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
//44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
// This is the same as the thing about it, but it includes a scaling factor
void cpu_iffty(float *real_image, float *imag_image, int size_x, int size_y, float*sin1, float*sin3, float *cos1,
		float *cos3) {
// Create some space for storing temporary values
	float *realOutBuffer = new float[size_y];
	float *imagOutBuffer = new float[size_y];
	float fft_real;
	float fft_imag;
	unsigned int n, m, x, y;
	float term;
	float *real_image_cp = new float[size_y];
	float *imag_image_cp = new float[size_y];

	float reduction_RealOutBuffer;
	float reduction_ImagOutBuffer;
	int eightX = size_x / 8;
	int eight7X = size_x - eightX;

	for (y = 0; y < size_y; y++) {

		for (m = 0; m < size_y; m++) {
			real_image_cp[m] = real_image[m * size_x + y];
			imag_image_cp[m] = imag_image[m * size_x + y];
		}

#pragma omp parallel for shared(real_image, imag_image, size_x, size_y, y, \
		realOutBuffer, imagOutBuffer, real_image_cp, imag_image_cp) \
        private(x, n, term, fft_real, fft_imag,reduction_RealOutBuffer, \
        reduction_ImagOutBuffer) schedule(dynamic, 4)

		for (x = 0; x < size_x; x++) {
			reduction_RealOutBuffer = 0.0f;
			reduction_ImagOutBuffer = 0.0f;

			for (n = 0; n < size_y; n++) {
				if (!(n >= eightX && n < eight7X)) {
					// Note that the negative sign goes away for the term
					term = 2 * PI * x * n / size_x;
//				fft_real = cos(term);
//				fft_imag = sin(term);

					fft_real = cos3[x * size_x + n];
					fft_imag = -sin3[x * size_x + n];

					reduction_RealOutBuffer += (real_image_cp[n] * fft_real) - (imag_image_cp[n] * fft_imag);
					reduction_ImagOutBuffer += (imag_image_cp[n] * fft_real) + (real_image_cp[n] * fft_imag);
				}
			}
			// Incorporate the scaling factor here
			realOutBuffer[x] = reduction_RealOutBuffer / size_x;
			imagOutBuffer[x] = reduction_ImagOutBuffer / size_x;
		}
// Write the buffer back to were the original values were
		for (int xx = 0; xx < size_x; xx++) {
			real_image[xx * size_x + y] = realOutBuffer[xx];
			imag_image[xx * size_x + y] = imagOutBuffer[xx];
		}
	}
// Reclaim some memory
	delete[] realOutBuffer;
	delete[] imagOutBuffer;
	delete[] real_image_cp;
	delete[] imag_image_cp;
}

void cpu_filter(float *real_image, float *imag_image, int size_x, int size_y) {
	int eightX = size_x / 8;
	int eight7X = size_x - eightX;
	int eightY = size_y / 8;
	int eight7Y = size_y - eightY;
	unsigned int x;
	unsigned int y;

#pragma omp parallel for shared(real_image, imag_image, size_x, size_y, eightX, eight7X, eightY, eight7Y) private(x, y)
	for (x = 0; x < size_x; x++) {
		for (y = 0; y < size_y; y++) {
			if (!(x < eightX && y < eightY) && !(x < eightX && y >= eight7Y) && !(x >= eight7X && y < eightY)
					&& !(x >= eight7X && y >= eight7Y)) {
				// Zero out these values
				real_image[y * size_x + x] = 0;
				imag_image[y * size_x + x] = 0;
			}
		}
	}
}

void cpu_filter_1(float *real_image, float *imag_image, int size_x, int size_y) {
	int eightX = size_x / 8;
	int eight7X = size_x - eightX;
	int eightY = size_y / 8;
	int eight7Y = size_y - eightY;
	unsigned int x;
	unsigned int y;

//#pragma omp parallel for shared(real_image, imag_image, size_x, size_y) private(x, y, eightX, eight7X, eightY, eight7Y)

#pragma omp parallel shared(real_image, imag_image, eightX, eightY, \
  							eight7X, size_x, eight7Y, size_y ) private(x,y)
#pragma omp sections nowait
	{
#pragma omp section
		{
			for (x = eightX; x < eight7X; x++) {
				for (y = 0; y < size_y; y++) {
					real_image[y * size_x + x] = 0;
					imag_image[y * size_x + x] = 0;
				}
			}
		}

#pragma omp section
		{
			for (x = 0; x < eightX; x++) {
				for (y = eightY; y < eight7Y; y++) {
					real_image[y * size_x + x] = 0;
					imag_image[y * size_x + x] = 0;
				}
			}
		}

#pragma omp section
		{
			for (x = eight7X; x < size_x; x++) {
				for (y = eightY; y < eight7Y; y++) {
					real_image[y * size_x + x] = 0;
					imag_image[y * size_x + x] = 0;
				}
			}
		}
	}
}

float imageCleaner(float *real_image, float *imag_image, int size_x, int size_y) {
// These are used for timing
	struct timeval tv1, tv2;
	struct timezone tz1, tz2;

// Start timing
	gettimeofday(&tv1, &tz1);

	float *sin1 = new float[size_y * size_y];
	float *sin3 = new float[size_y * size_y];

	float *cos1 = new float[size_x * size_y];
	float *cos3 = new float[size_x * size_y];
	precalculate(sin1, sin3, cos1, cos3, size_x, size_y);
// Perform fft with respect to the x direction
	cpu_fftx(real_image, imag_image, size_x, size_y, sin1, sin3, cos1, cos3);
// Perform fft with respect to the y direction
	cpu_ffty(real_image, imag_image, size_x, size_y, sin1, sin3, cos1, cos3);

// Filter the transformed image
//	cpu_filter(real_image, imag_image, size_x, size_y);

// Perform an inverse fft with respect to the x direction
	cpu_ifftx(real_image, imag_image, size_x, size_y, sin1, sin3, cos1, cos3);
// Perform an inverse fft with respect to the y direction
	cpu_iffty(real_image, imag_image, size_x, size_y, sin1, sin3, cos1, cos3);

// End timing
	gettimeofday(&tv2, &tz2);

// Compute the time difference in micro-seconds
	float execution = ((tv2.tv_sec - tv1.tv_sec) * 1000000 + (tv2.tv_usec - tv1.tv_usec));
// Convert to milli-seconds
	execution /= 1000;
// Print some output
	printf("OPTIMIZED IMPLEMENTATION STATISTICS:\n");
	printf("  Optimized Kernel Execution Time: %f ms\n\n", execution);
	return execution;
}

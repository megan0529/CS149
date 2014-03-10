We found that the discrete fourier transform of a picture is in fact the result of matrix multiplication: matrix A being the original picture, matrix B being the sin and cos matrix. So we implmented a kernel of matrix multiplication in cuda. We took advantage of shared memory by blocking the matrix. Thus the amount of memory access is equal to the (matrix width / block width).

One thing to be noticed is that in the FT_Kernel the result has to be write to another place, instead of write back to the orignal place. This is due to the fact that other thread blocks might need this matrix data, so you can't overwirte it.

Another trick is when doing the column wise fourier transform, it needs to do transpose of the input matrix. This is done by storing the result of row-wise transform in transposed order. Thus no explicit transpose is needed. 

The final speed up for the 1024x1024 input is:
	CUDA IMPLEMENTATION STATISTICS:
	Host to Device Transfer Time: 2.787520 ms
	Kernel(s) Execution Time: 111.720673 ms
	Device to Host Transfer Time: 2.698048 ms
	Total CUDA Execution Time: 117.206245 ms

TOTAL SPEEDUP: 4136.980957


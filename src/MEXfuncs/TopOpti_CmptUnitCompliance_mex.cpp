#include "mex.h"
#include <omp.h>
#include <stdlib.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" TopOpti_CmptUnitCompliance_mex.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "Four inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "One output required.");
    }

    // Retrieve the inputs
    double *U = mxGetPr(prhs[0]);
    int32_T *eNodMat = (int32_T *)mxGetData(prhs[1]);
    double *Ke = mxGetPr(prhs[2]);
    double *blockSizeInput = mxGetPr(prhs[3]);

    // Validate the block size input
    mwSize blockSize = (mwSize)blockSizeInput[0];
    if (blockSize <= 0) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidBlockSize", "Block size must be a positive integer.");
    }

    // Validate input dimensions
    if (!mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input U must be a column vector.");
    }
    if (!mxIsInt32(prhs[1]) || mxGetN(prhs[1]) != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input eNodMat must be an Mx8 int32 matrix.");
    }
    if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 24 || mxGetN(prhs[2]) != 24) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNot24x24", "Input Ke must be a 24x24 double matrix.");
    }

    // Retrieve dimensions
    mwSize numDOFs = mxGetM(prhs[0]);
    mwSize numElements = mxGetM(prhs[1]);

    // Create output array Y with dimensions [numElements, 1]
    plhs[0] = mxCreateDoubleMatrix(numElements, 1, mxREAL);
    double *Y = mxGetPr(plhs[0]);
	int i;
    // Initialize Y to zeros
    #pragma omp parallel for schedule(static)
    for (i = 0; i < numElements; ++i) {
        Y[i] = 0.0;
    }

    // Process eNodMat in blocks
    for (mwSize block_start = 0; block_start < numElements; block_start += blockSize) {
        mwSize block_end = block_start + blockSize;
        if (block_end > numElements) {
            block_end = numElements;
        }

        mwSize block_size = block_end - block_start;

        // Create Umat matrix with dimensions [block_size, 24]
        double *Umat = (double *)mxMalloc(block_size * 24 * sizeof(double));

        // Populate Umat using eNodMat and U
        #pragma omp parallel for schedule(static)
        for (i = 0; i < block_size; ++i) {
            mwSize global_index = block_start + i;
            for (int j = 0; j < 8; ++j) {
                int baseIndex = 3 * (eNodMat[global_index + j * numElements] - 1);
                for (int k = 0; k < 3; ++k) {
                    Umat[i * 24 + 3 * j + k] = U[baseIndex + k];
                }
            }
        }

        // Create X matrix with dimensions [block_size, 24]
        double *X = (double *)mxMalloc(block_size * 24 * sizeof(double));

        // Multiply each row of Umat with Ke and scale by the corresponding elements in E
        #pragma omp parallel for schedule(static)
        for (i = 0; i < block_size; ++i) {
            for (int j = 0; j < 24; ++j) {
				double val = 0;
                for (int k = 0; k < 24; ++k) {
					val += Umat[i * 24 + k] * Ke[k * 24 + j];
                }
				X[i * 24 + j] = val;
            }
        }
		
		// Multiply Umat * X
		#pragma omp parallel for schedule(static)
		for (i = 0; i < block_size; ++i) {
			mwSize global_index = block_start + i;
			double val = 0;
			for (int j = 0; j < 24; ++j) {
				val += Umat[i * 24 + j]*X[i * 24 + j];
			}
			Y[global_index] = val;
		}

        // Free temporary matrices
        mxFree(Umat);
        mxFree(X);
    }
}

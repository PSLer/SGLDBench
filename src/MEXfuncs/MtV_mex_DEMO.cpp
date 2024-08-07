#include "mex.h"
#include <omp.h>
//mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" MtV_mex_DEMO.cpp
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "Four inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "One output required.");
    }

    // Input validation
    if (!mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input U must be a column vector.");
    }
    if (!mxIsInt32(prhs[1]) || mxGetN(prhs[1]) != 24) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input eDofMat must be an Mx24 int32 matrix.");
    }
    if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 24 || mxGetN(prhs[2]) != 24) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNot24x24", "Input Ke must be a 24x24 double matrix.");
    }
    if (!mxIsDouble(prhs[3]) || mxGetN(prhs[3]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input E must be a column vector.");
    }

    // Retrieve the inputs
    double *U = mxGetPr(prhs[0]);
    int32_T *eDofMat = (int32_T *)mxGetData(prhs[1]);
    double *Ke = mxGetPr(prhs[2]);
    double *E = mxGetPr(prhs[3]);

    mwSize numDOFs = mxGetM(prhs[0]);
    mwSize numElements = mxGetM(prhs[1]);

    // Create output array Y with dimensions [numDOFs, 1]
    plhs[0] = mxCreateDoubleMatrix(numDOFs, 1, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    // Initialize Y to zeros
	int i;
    #pragma omp parallel for
    for (i = 0; i < numDOFs; ++i) {
        Y[i] = 0.0;
    }

    // Create temporary matrix Umat with dimensions [numElements, 24]
    double *Umat = (double *)mxMalloc(numElements * 24 * sizeof(double));

    // Indexing the elements in U into Umat
    #pragma omp parallel for
    for (i = 0; i < numElements; ++i) {
        for (int j = 0; j < 24; ++j) {
            Umat[i * 24 + j] = U[eDofMat[i + j * numElements] - 1]; // 1-based to 0-based index
        }
    }

    // Create temporary matrix X with dimensions [numElements, 24]
    double *X = (double *)mxMalloc(numElements * 24 * sizeof(double));

    // Multiply each row of Umat with Ke and scale by the corresponding elements in E
    #pragma omp parallel for
    for (i = 0; i < numElements; ++i) {
        for (int j = 0; j < 24; ++j) {
            X[i * 24 + j] = 0.0;
            for (int k = 0; k < 24; ++k) {
                X[i * 24 + j] += Umat[i * 24 + k] * Ke[k * 24 + j];
            }
            X[i * 24 + j] *= E[i];
        }
    }

    // Accumulate the elements in X into Y
    #pragma omp parallel for
    for ( i = 0; i < numElements; ++i) {
        for (int j = 0; j < 24; ++j) {
            #pragma omp atomic
            Y[eDofMat[i + j*numElements] - 1] += X[i * 24 + j]; // 1-based to 0-based index
        }
    }

    // Free temporary matrices
    mxFree(Umat);
    mxFree(X);
}

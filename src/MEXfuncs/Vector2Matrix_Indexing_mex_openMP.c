#include "mex.h"
#include "omp.h"
//Compelling
//mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Vector2Matrix_Indexing_mex_openMP.c
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for proper number of arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:indexing_mex:nrhs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:indexing_mex:nlhs", "One output required.");
    }
    
    // Ensure mapVec is of int32 type
    if (!mxIsInt32(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:indexing_mex:notInt32", "mapVec must be of type int32.");
    }
	
    // Input arguments
    double *V = mxGetPr(prhs[0]);
    int32_T *mapVec = (int32_T *)mxGetData(prhs[1]);
    
    // Dimensions of mapVec
    mwSize rows = mxGetM(prhs[1]);
    mwSize cols = mxGetN(prhs[1]);
    
    // Create the output matrix
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double *M = mxGetPr(plhs[0]);
    
    // Perform the mapping using OpenMP parallel for loop
	/*mwSize j, i;
    #pragma omp parallel for collapse(2) private(i)
    for (j = 0; j < cols; ++j) {
        for (i = 0; i < rows; ++i) {
            mwSize idx = (mwSize)mapVec[i + j * rows] - 1; // MATLAB uses 1-based indexing
            M[i + j * rows] = V[idx];
        }
    }*/
	int jx;
    #pragma omp parallel for
    for (jx = 0; jx < cols; ++jx) {
		mwSize j = (mwSize)jx;
        for (mwSize i = 0; i < rows; ++i) {
            mwSize idx = (mwSize)mapVec[i + j * rows] - 1; // MATLAB uses 1-based indexing
            M[i + j * rows] = V[idx];
        }
    }
}	
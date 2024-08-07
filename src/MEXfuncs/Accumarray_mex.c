#include "mex.h"
#include <stddef.h> // for NULL

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for proper number of arguments
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:accumulation_mex:nrhs", "Three inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:accumulation_mex:nlhs", "One output required.");
    }
    
    // Retrieve input arguments
    int32_T *A = (int32_T *)mxGetData(prhs[0]);
    double *B = mxGetPr(prhs[1]);
    mwSize M = mxGetScalar(prhs[2]); // Extract M directly from [M, 1]
    
    mwSize n = mxGetNumberOfElements(prhs[0]);
    
    // Ensure A and B are column vectors
    if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetNumberOfDimensions(prhs[1]) != 2) {
        mexErrMsgIdAndTxt("MyToolbox:accumulation_mex:invalidInput", "Inputs A and B must be column vectors.");
    }
    if (mxGetM(prhs[0]) != n || mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:accumulation_mex:invalidInput", "Input A must be a column vector.");
    }
    if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:accumulation_mex:invalidInput", "Input B must be a column vector.");
    }

    // Create the output matrix C with dimension [M, 1]
    plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);
    double *C = mxGetPr(plhs[0]);
    
    // Initialize C with zeros
    for (mwSize i = 0; i < M; ++i) {
        C[i] = 0.0;
    }
    
    // Perform the accumulation
    for (mwSize i = 0; i < n; ++i) {
        mwIndex idx = A[i] - 1; // Convert 1-based index to 0-based
        if (idx < 0 || idx >= M) {
            // Handle out-of-bounds index
            mexErrMsgIdAndTxt("MyToolbox:accumulation_mex:indexOutOfBounds", "Index %d (0-based: %d) is out of bounds. Valid indices are 1 to %zu.", A[i], idx, M);
        }
        C[idx] += B[i];
    }
}

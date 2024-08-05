#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:indexing_mex:nrhs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:indexing_mex:nlhs", "One output required.");
    }

    if (!mxIsInt32(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:indexing_mex:notInt32", "mapVec must be of type int32.");
    }

    double *V = mxGetPr(prhs[0]);
    int32_T *mapVec = (int32_T *)mxGetData(prhs[1]);
    
    mwSize rows = mxGetM(prhs[1]);
    mwSize cols = mxGetN(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double *M = mxGetPr(plhs[0]);
    
    for (mwSize i = 0; i < rows; ++i) {
        for (mwSize j = 0; j < cols; ++j) {
            mwSize idx = (mwSize)mapVec[i + j * rows] - 1;
            M[i + j * rows] = V[idx];
        }
    }
}
// mexFunction.cpp
#include "mex.h"
#include <omp.h>
#include "MMAseq.h"
// mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" MMA_mex.cpp MMAseq.cc
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Input validation (optional)
    if (nrhs != 10) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "10 inputs required.");
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "3 output required.");
    }
	int numConstraints = (int) mxGetScalar(prhs[0]); //m
	int numVariables = (int) mxGetScalar(prhs[1]); //n
	double *xVal = mxGetPr(prhs[2]);
	double *xmin = mxGetPr(prhs[3]);
	double *xmax = mxGetPr(prhs[4]);	
	double *xold1 = mxGetPr(prhs[5]);
	double *xold2 = mxGetPr(prhs[6]);		
	double *df0dx = mxGetPr(prhs[7]);
	double *fval = mxGetPr(prhs[8]);
	double *dfdx = mxGetPr(prhs[9]);
	
	plhs[0] = mxCreateDoubleMatrix(numVariables, 1, mxREAL);
	double *xnew = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(numVariables, 1, mxREAL);
	double *xold1_out = mxGetPr(plhs[1]);
	plhs[2] = mxCreateDoubleMatrix(numVariables, 1, mxREAL);
	double *xold2_out = mxGetPr(plhs[2]);
	
	int i;
	double *xVal_tmp = (double*) mxMalloc(numVariables * sizeof(double));
	MMAseq *mma = new MMAseq(numVariables, numConstraints, 0.0, 100.0, 0.0);
	
	#pragma omp parallel for
	for (i = 0; i < numVariables; ++i) {
		xVal_tmp[i] = xVal[i];
		mma->xo1[i] = xold1[i];
		mma->xo2[i] = xold2[i];
	}

	//mma->Update(&xval_[0], &df0dx_[0], &fval_[0], &dfdx_[0], &xmin_[0], &xmax_[0]);
	mma->Update(xVal_tmp, df0dx, fval, dfdx, xmin, xmax);
	#pragma omp parallel for
	for (i = 0; i < numVariables; ++i) {
		xnew[i] = xVal_tmp[i];
		xold1_out[i] = mma->xo1[i];
		xold2_out[i] = mma->xo2[i];
	}

	mxFree(xVal_tmp);
	delete mma;
}

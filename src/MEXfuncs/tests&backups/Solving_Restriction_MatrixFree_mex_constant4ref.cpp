#include "mex.h"
#include <omp.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_Restriction_MatrixFree_mex_constant4ref.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "Five inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "One output required.");
    }

    // Retrieve the inputs
    double *rFiner = mxGetPr(prhs[0]);
    int32_T *transferMat = (int32_T *)mxGetData(prhs[1]); //125 or 27 * numElements
    double *multiGridOperatorRI = mxGetPr(prhs[2]);
    int32_T *eNodMat = (int32_T *)mxGetData(prhs[3]);
    mwSize numNodes = (mwSize)mxGetScalar(prhs[4]);

    // Retrieve dimensions
    mwSize nFiner = mxGetM(prhs[0]); // Number of rowsTransferMat in rFiner
    
    mwSize anotherDimTransferMat = mxGetM(prhs[1]); // Number of columns in transferMat (27 or 125)
    mwSize numMultiGridRows = mxGetM(prhs[2]); // Number of rowsTransferMat in multiGridOperatorRI (125 or 27)
    mwSize numMultiGridCols = mxGetN(prhs[2]); // Number of columns in multiGridOperatorRI (8)
	mwSize numElements = mxGetN(prhs[1]); // Number of rowsTransferMat in transferMat
    
	mwSize rowsTransferMat = mxGetM(prhs[1]); // number of elements 
    /* mwSize colsTransferMat = mxGetN(prhs[1]); 	 */
    // Validate input dimensions
    if (!mxIsInt32(prhs[1]) || !mxIsInt32(prhs[3])) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotInt32", "transferMat and eNodMat must be of type int32.");
    }
    if (anotherDimTransferMat != numMultiGridRows || numMultiGridCols != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidMultiGridOperatorRI", "multiGridOperatorRI must have 8 columns.");
    }
	int j;
if (1) //Tests 
{
	plhs[0] = mxCreateDoubleMatrix(numNodes, 1, mxREAL);
	double *Y = mxGetPr(plhs[0]);
	
	double *Umat = (double *)mxMalloc(numElements * anotherDimTransferMat * sizeof(double));

	// Populate Umat using transferMat and rFiner
    #pragma omp parallel for
    for (j = 0; j < numElements; ++j) {
        for (mwSize i = 0; i < anotherDimTransferMat; ++i) {
			Umat[i + j * anotherDimTransferMat] = rFiner[transferMat[i + j * anotherDimTransferMat] - 1];
        }
    }

    // Create X matrix with dimensions [numElements, 24]
    double *X = (double *)mxMalloc(numElements * 8 * sizeof(double));
		
	// Multiply each column of Umat with multiGridOperatorRI
	#pragma omp parallel for schedule(static)
	for (j = 0; j < numElements; ++j) {
		for (mwSize i = 0; i < 8; ++i) {
			X[i * numElements + j] = 0.0;
				for (int k = 0; k < anotherDimTransferMat; ++k) {
					double valueAnotherDimTransferMat =  multiGridOperatorRI[i*anotherDimTransferMat + k];
					if (0!=valueAnotherDimTransferMat) {
						X[i * numElements + j] += Umat[k + j * anotherDimTransferMat] * valueAnotherDimTransferMat;
					}
					
				
/* 					//Previous
					X[i * numElements + j] += Umat[k + j * anotherDimTransferMat] * multiGridOperatorRI[i*anotherDimTransferMat + k];
					 */
				}
		}
	}
	
	// Accumulate the elements in X into Y through eNodMat
	#pragma omp parallel for schedule(static)
	for (j = 0; j < numElements; ++j) {
		for (mwSize i = 0; i < 8; ++i) {
			int baseIndex = eNodMat[j + i * numElements] - 1;
			#pragma omp atomic
			Y[baseIndex] += X[i * numElements + j];
		}
	}
	
    // Free temporary matrices
    mxFree(Umat);
	mxFree(X);	
}
else 
{
	
	//Test 2
	plhs[0] = mxCreateDoubleMatrix(numElements, 8, mxREAL);
	double *X = mxGetPr(plhs[0]);
	
	double *Umat = (double *)mxMalloc(numElements * anotherDimTransferMat * sizeof(double));

	// Populate Umat using transferMat and rFiner
    #pragma omp parallel for
    for (j = 0; j < numElements; ++j) {
        for (mwSize i = 0; i < anotherDimTransferMat; ++i) {
			Umat[i + j * anotherDimTransferMat] = rFiner[transferMat[i + j * anotherDimTransferMat] - 1];
        }
    }
	
	// Multiply each column of Umat with multiGridOperatorRI
	#pragma omp parallel for schedule(static)
	for (j = 0; j < numElements; ++j) {
		for (mwSize i = 0; i < 8; ++i) {
			X[i * numElements + j] = 0.0;
				for (int k = 0; k < anotherDimTransferMat; ++k) {
					X[i * numElements + j] += Umat[k + j * anotherDimTransferMat] * multiGridOperatorRI[i*anotherDimTransferMat + k];
					
				}
		}
	}
	
    // Free temporary matrices
    mxFree(Umat);	
}

}

/* 	//Test 1
	plhs[0] = mxCreateDoubleMatrix(anotherDimTransferMat, numElements, mxREAL);
	double *Umat = mxGetPr(plhs[0]);

    #pragma omp parallel for
    for (j = 0; j < numElements; ++j) {
        for (mwSize i = 0; i < anotherDimTransferMat; ++i) {
			Umat[i + j * anotherDimTransferMat] = rFiner[transferMat[i + j * anotherDimTransferMat] - 1];
        }
    } */
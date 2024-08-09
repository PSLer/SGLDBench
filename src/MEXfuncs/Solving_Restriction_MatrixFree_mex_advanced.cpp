#include "mex.h"
#include <stdint.h>
#include <omp.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_Restriction_MatrixFree_mex_advanced.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "Five inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "One output required.");
    }

    // Retrieve the inputs
    double *rFiner = mxGetPr(prhs[0]);
    int32_T *transferMat = (int32_T *)mxGetData(prhs[1]); //125 or 27 * numElements
    /* double *multiGridOperatorRI = mxGetPr(prhs[2]); */
	const mxArray *multiGridOperatorRI = prhs[2]; // Sparse matrix
    int32_T *eNodMat = (int32_T *)mxGetData(prhs[3]);
    mwSize numNodes = (mwSize)mxGetScalar(prhs[4]);
	mwSize intermediateNumNodes = (mwSize)mxGetScalar(prhs[5]);
	int32_T *solidNodeMapCoarser2Finer = (int32_T *)mxGetData(prhs[6]);
	double *transferMatCoeffi = mxGetPr(prhs[7]);
	
    // Retrieve dimensions
    mwSize nFiner = mxGetM(prhs[0]); 
	mwSize colFiner = mxGetN(prhs[0]); //constant 3
    
    mwSize anotherDimTransferMat = mxGetM(prhs[1]); // Number of columns in transferMat (27 or 125)
    mwSize numMultiGridRows = mxGetM(prhs[2]); // Number of rowsTransferMat in multiGridOperatorRI (125 or 27)
    mwSize numMultiGridCols = mxGetN(prhs[2]); // Number of columns in multiGridOperatorRI (8)
	mwSize numElements = mxGetN(prhs[1]); // Number of rowsTransferMat in transferMat
	
	
    // Validate input dimensions
    if (!mxIsInt32(prhs[1]) || !mxIsInt32(prhs[3])) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotInt32", "transferMat and eNodMat must be of type int32.");
    }
    if (anotherDimTransferMat != numMultiGridRows || numMultiGridCols != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidMultiGridOperatorRI", "multiGridOperatorRI must have 8 columns.");
    }
    if (mxGetN(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:solidNodeMapCoarser2Finer", "solidNodeMapCoarser2Finer must have 1 columns.");
    }
	int j;
	plhs[0] = mxCreateDoubleMatrix(numNodes, colFiner, mxREAL);
	double *Y = mxGetPr(plhs[0]);
    
	// Retrieve sparse matrix data
    mwIndex *ir = mxGetIr(multiGridOperatorRI);
    mwIndex *jc = mxGetJc(multiGridOperatorRI);
    const double *sparseData = mxGetPr(multiGridOperatorRI);
	
	double *Umat = (double *)mxMalloc(numElements * anotherDimTransferMat * sizeof(double));

    // Create X matrix with dimensions [numElements, 24]
    double *X = (double *)mxMalloc(numElements * 8 * sizeof(double));
	double *rFiner1Child = (double *)mxMalloc(intermediateNumNodes * sizeof(double));
	double *YChild = (double *)mxMalloc(numNodes * sizeof(double));
	
	for (int ss = 0; ss < colFiner; ++ss) {
		//Initialize child variables
		#pragma omp parallel for
		for (j = 0; j < intermediateNumNodes; ++j) {
			rFiner1Child[j] = 0.0;
		}

		#pragma omp parallel for
		for (j = 0; j < nFiner; ++j) {
			rFiner1Child[solidNodeMapCoarser2Finer[j]-1] = rFiner[ss*nFiner + j];
		}
		
		#pragma omp parallel for
		for (j = 0; j < intermediateNumNodes; ++j) {
			rFiner1Child[j] /= transferMatCoeffi[j];
		}		

		#pragma omp parallel for
		for (j = 0; j < numNodes; ++j) {
			YChild[j] = 0.0;
		}
		
		// Populate Umat using transferMat and rFiner1Child
		#pragma omp parallel for
		for (j = 0; j < numElements; ++j) {
			for (mwSize i = 0; i < anotherDimTransferMat; ++i) {
				Umat[i + j * anotherDimTransferMat] = rFiner1Child[transferMat[i + j * anotherDimTransferMat] - 1];
			}
		}
		
		// Multiply each column of Umat with multiGridOperatorRI
		#pragma omp parallel for
		for ( j = 0; j < 8; ++j) {
			for (mwSize i = 0; i < numElements; ++i) {
				double sum = 0.0;
				// Iterate over non-zero entries in column i of multiGridOperatorRI
				for (mwIndex k = jc[j]; k < jc[j + 1]; ++k) {
					mwIndex row = ir[k];
					sum += Umat[row + i * anotherDimTransferMat] * sparseData[k];
				}
				X[i + j * numElements] = sum;
			}
		}		

		// Accumulate the elements in X into Y through eNodMat
		#pragma omp parallel for schedule(static)
		for (j = 0; j < numElements; ++j) {
			for (mwSize i = 0; i < 8; ++i) {
				int baseIndex = eNodMat[j + i * numElements] - 1;
				#pragma omp atomic
				YChild[baseIndex] += X[i * numElements + j];
			}
		}

		// store to Y		
		#pragma omp parallel for
		for (j = 0; j < numNodes; ++j) {
			Y[j+ss*numNodes] = YChild[j];
		}
	}
	
    // Free temporary matrices
    mxFree(Umat);
	mxFree(X);	
	mxFree(YChild);
	mxFree(rFiner1Child);
}

#include "mex.h"
#include <stdint.h>
#include <omp.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_Interpolation_MatrixFree_mex.cpp"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "Five inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "One output required.");
    }

    // Retrieve the inputs
    double *xCoarser = mxGetPr(prhs[0]);
    int32_T *transferMat = (int32_T *)mxGetData(prhs[1]); //125 or 27 * numElements
    /* double *multiGridOperatorRI = mxGetPr(prhs[2]); */
	const mxArray *multiGridOperatorRI = prhs[2]; // Sparse matrix
    int32_T *eNodMat = (int32_T *)mxGetData(prhs[3]);
    mwSize numNodes = (mwSize)mxGetScalar(prhs[4]);
	mwSize intermediateNumNodes = (mwSize)mxGetScalar(prhs[5]);
	int32_T *solidNodeMapCoarser2Finer = (int32_T *)mxGetData(prhs[6]);
	double *transferMatCoeffi = mxGetPr(prhs[7]);
	
    // Retrieve dimensions
    mwSize nCoarser = mxGetM(prhs[0]); 
	mwSize colFiner = mxGetN(prhs[0]); //constant 3
    
    mwSize anotherDimTransferMat = mxGetM(prhs[1]); // Number of columns in transferMat (27 or 125)
    mwSize numMultiGridRows = mxGetM(prhs[2]); // Number of rowsTransferMat in multiGridOperatorRI (125 or 27)
    mwSize numMultiGridCols = mxGetN(prhs[2]); // Number of columns in multiGridOperatorRI (8)
	mwSize numElements = mxGetM(prhs[3]); // Number of rowsTransferMat in transferMat
	
	
    // Validate input dimensions
    if (!mxIsInt32(prhs[1]) || !mxIsInt32(prhs[3])) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotInt32", "transferMat and eNodMat must be of type int32.");
    }
/*     if (anotherDimTransferMat != numMultiGridCols || numMultiGridRows != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidMultiGridOperatorRI", "multiGridOperatorRI must have 8 columns.");
    } */
    if (mxGetN(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:solidNodeMapCoarser2Finer", "solidNodeMapCoarser2Finer must have 1 columns.");
    }
	int j;

	//default output
	plhs[0] = mxCreateDoubleMatrix(numNodes, colFiner, mxREAL);
	double *Y = mxGetPr(plhs[0]);	
		
	// Retrieve sparse matrix data
    mwIndex *ir = mxGetIr(multiGridOperatorRI);
    mwIndex *jc = mxGetJc(multiGridOperatorRI);
	mwSize numRowsMultiGridOperatorRI = mxGetM(prhs[1]); // Should be 125
    mwSize numColsMultiGridOperatorRI = mxGetN(prhs[1]); // Should be 8
    const double *sparseData = mxGetPr(multiGridOperatorRI);
	
	double *Umat = (double *)mxMalloc(numElements * 8 * sizeof(double));

    // Create X matrix with dimensions [numElements, 24]
    double *X = (double *)mxMalloc(numElements * anotherDimTransferMat * sizeof(double));
	double *rCoarseChild = (double *)mxMalloc(nCoarser * sizeof(double));
	double *YChild = (double *)mxMalloc(intermediateNumNodes * sizeof(double));
	double *YChildCompact = (double *)mxMalloc(numNodes * sizeof(double));
	
	for (int ss = 0; ss < colFiner; ++ss) {
		//Initialize child variables

		#pragma omp parallel for
		for (j = 0; j < nCoarser; ++j) {
			rCoarseChild[j] = xCoarser[ss*nCoarser + j];
		}

		#pragma omp parallel for
		for (j = 0; j < intermediateNumNodes; ++j) {
			YChild[j] = 0.0;
		}
		
		// Populate Umat using eNodMat and rCoarseChild
 		#pragma omp parallel for
		for (j = 0; j < numElements; ++j) {
			for (mwSize i = 0; i < 8; ++i) {
				Umat[i + j * 8] = rCoarseChild[eNodMat[j + i * numElements] - 1];
			}
		}

		#pragma omp parallel for
		for (j = 0; j < anotherDimTransferMat * numElements; j++) {
			X[j] = 0.0;
		}
		
		// Perform multiplication X = multiGridOperatorRI * Umat
		#pragma omp parallel for
		for (j = 0; j < numElements; j++) {  // Loop over columns of Umat
			for (mwSize k = 0; k < 8; k++) {  // Loop over columns of multiGridOperatorRI
				for (mwSize idx = jc[k]; idx < jc[k + 1]; idx++) {  // Loop over non-zero entries in column k
					mwSize i = ir[idx];  // Row index of non-zero element
					X[i + j * anotherDimTransferMat] += sparseData[idx] * Umat[k + j * 8];
				}
			}
		}	

		// Accumulate the elements in X into Y through eNodMat
		#pragma omp parallel for schedule(static)
		for (j = 0; j < numElements; ++j) {
			for (mwSize i = 0; i < anotherDimTransferMat; ++i) {
				int baseIndex = transferMat[i + j * anotherDimTransferMat] - 1;
				#pragma omp atomic
				YChild[baseIndex] += X[i + j * anotherDimTransferMat];
			}
		}
		
		#pragma omp parallel for
		for (j = 0; j < intermediateNumNodes; ++j) {
			YChild[j] /= transferMatCoeffi[j];
		}

		//Extract the compact YChild
		#pragma omp parallel for
		for (j = 0; j < numNodes; ++j) {
			YChildCompact[j] = YChild[solidNodeMapCoarser2Finer[j]-1];
		}

		// store to Y		
		#pragma omp parallel for
		for (j = 0; j < numNodes; ++j) {
			Y[j+ss*numNodes] = YChildCompact[j];
		}
	}
	
    // Free temporary matrices
    mxFree(Umat);
	mxFree(X);	
	mxFree(YChild);
	mxFree(YChildCompact);
}

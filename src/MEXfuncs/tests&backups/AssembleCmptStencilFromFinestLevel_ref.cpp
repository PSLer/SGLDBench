#include "mex.h"
#include "matrix.h"
#include <omp.h>

/* Declaration of the independent functions */
void AccumulationSmall(double *sK, double *tmpK, int *localMapping, mwSize nPerParentEleDOFs, mwSize perParentEleSons);
void ComputeIntermediateResult(double *B, const mwIndex *rowIndices, const mwIndex *colPointers, 
                               const double *values, const double *tmpK, mwSize nPerParentEleDOFs, mwSize nRows);
void ComputeFinalResult(double *iKeCoarser, const double *B, const mwIndex *rowIndices, 
                        const mwIndex *colPointers, const double *values, mwSize nPerParentEleDOFs, mwSize nRows);

/* Main MEX function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Variable declarations */
    double *iKe;                   /* 576x1 matrix */
    double *eleModulus;           /* Row vector with nElesFinest elements */
    int *elementUpwardMap;        /* nElesCurrent x perParentEleSons matrix (int32) */
    const mxArray *interpolatingKe; /* Sparse matrix with dimensions nPerParentEleDOFs x 24 */
    int *localMapping;            /* Column vector with 24*24*perParentEleSons elements (int32) */
    mwSize numProjectNodes;          /* Scalar int32 value */
    
    /* Get input arguments */
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Six inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "One output required.");
    }

    iKe = mxGetPr(prhs[0]);                    /* Input column vector iKe (576x1) */
    eleModulus = mxGetPr(prhs[1]);            /* Row vector eleModulus (nElesFinest) */
    elementUpwardMap = (int *) mxGetData(prhs[2]); /* Matrix elementUpwardMap (nElesCurrent x perParentEleSons) as int32 */
    interpolatingKe = prhs[3];         /* Sparse matrix interpolatingKe (nPerParentEleDOFs x 24) */
    localMapping = (int *) mxGetData(prhs[4]); /* Column vector localMapping (24*24*perParentEleSons x 1) as int32 */
    numProjectNodes = (int) mxGetScalar(prhs[5]);  /* Scalar numProjectNodes as int32 */

    /* Compute derived dimensions */
    mwSize nPerParentEleDOFs = 3 * numProjectNodes;
    mwSize perParentEleSons = mxGetN(prhs[2]);  /* Number of columns in elementUpwardMap */
    mwSize nElesCurrent = mxGetM(prhs[2]);  /* Number of rows in elementUpwardMap */

    /* Initialize Ks as a full zero array */
    mwSize dims[3] = {24, 24, nElesCurrent};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *Ks = mxGetPr(plhs[0]);

    /* Declare and allocate sK outside of the loops */
    mwSize sK_size = 576 * perParentEleSons;
    double *sK = (double*) mxMalloc(sK_size * sizeof(double)); 

    /* Declare tmpK */
	double *tmpK = (double*) mxMalloc(nPerParentEleDOFs * nPerParentEleDOFs * sizeof(double));
	
    /* Extract sparse matrix information */
    mwIndex *rowIndices = mxGetIr(interpolatingKe);  /* Row indices of non-zero elements */
    mwIndex *colPointers = mxGetJc(interpolatingKe); /* Column pointers for the sparse matrix */
    double *values = mxGetPr(interpolatingKe);       /* Non-zero values in the sparse matrix */
    /* Allocate memory for intermediate result B (nRows x nPerParentEleDOFs) */
    double *B = (double *)mxCalloc(24 * nPerParentEleDOFs, sizeof(double));
    /* Create an output matrix for iKeCoarser (24x24) */
	mwSize nRows = 24; 
    double *iKeCoarser = (double *)mxMalloc(nRows * nRows * sizeof(double));

    int i;
	//#pragma omp parallel for
    //#pragma omp parallel for shared(Ks, sK, iKe, tmpK, iKeCoarser, localMapping) /* Shared variables among threads */
	//#pragma omp parallel for schedule(static)
	for (i = 0; i < nElesCurrent; i++) {
        /* Zero-initialize sK before entering the inner loop */
        for (mwSize k = 0; k < sK_size; k++) {
            sK[k] = 0.0;
        }

        /* Inner loop traverses the columns of elementUpwardMap */
        for (mwSize j = 0; j < perParentEleSons; j++) {
            int idx = elementUpwardMap[i + j * nElesCurrent]; /* Access the elementUpwardMap value */

            /* Skip if idx is 0 (unrelated value) */
            if (idx != 0) {
				/* Access eleModulus using idx - 1 (MATLAB indexing starts from 1, C from 0) */
				double modulus_value = eleModulus[idx - 1];
	
				/* Multiply iKe by modulus_value and store it in the j-th column of sK */
				for (mwSize k = 0; k < 576; k++) {
					sK[k + j * 576] = iKe[k] * modulus_value;
				}                
            }
        }

        /* Accumulate the values in sK using localMapping and reshape to sparse format */
		for (int j = 0; j < nPerParentEleDOFs * nPerParentEleDOFs; j++) {
			tmpK[j] = 0.0;
		}		
        /* AccumulationSmall(sK, tmpK, localMapping, nPerParentEleDOFs, perParentEleSons); */
		/* Perform the accumulation */
		for (int j = 0; j < perParentEleSons; j++) {
			for (mwSize k = 0; k < 576; k++) {
				int idx = localMapping[k + j * 576] - 1;
				//#pragma omp atomic
				tmpK[idx] += sK[k + j * 576];
			}
		}
	
        /* Compute iKeCoarser_dense using sparse matrix operations */
		for (mwSize j = 0; j < nRows * nPerParentEleDOFs; j++) {
			B[j] = 0.0;
		}		
		//ComputeIntermediateResult(B, rowIndices, colPointers, values, tmpK, nPerParentEleDOFs, nRows);
		for (mwSize col = 0; col < nRows; col++) {
			for (mwIndex k = colPointers[col]; k < colPointers[col + 1]; k++) {
				mwSize ss = rowIndices[k];
				double val = values[k];
				for (mwSize row = 0; row < nPerParentEleDOFs; row++) {
					B[col * nPerParentEleDOFs + row] += val * tmpK[row * nPerParentEleDOFs + ss];
				}
			}
		}		
		for (mwSize j = 0; j < nRows * nRows; j++) {
			iKeCoarser[j] = 0.0;
		}		
		/* ComputeFinalResult(iKeCoarser, B, rowIndices, colPointers, values, nPerParentEleDOFs, nRows); */
		for (mwSize col = 0; col < nRows; col++) {
			for (mwIndex k = colPointers[col]; k < colPointers[col + 1]; k++) {
				mwSize ss = rowIndices[k];
				double val = values[k];
				for (mwSize row = 0; row < nRows; row++) {
					iKeCoarser[col * nRows + row] += B[row * nPerParentEleDOFs + ss] * val;
				}
			}
		}
	
        /* Copy the dense result into the i-th page of Ks */
        for (int row = 0; row < 24; row++) {
            for (int col = 0; col < 24; col++) {
                Ks[row + col * 24 + i * 576] = iKeCoarser[row + col * 24];
            }
        }

        /* Free memory for tmpK and iKeCoarser_dense if necessary */
/*         mxDestroyArray(tmpK);
        mxDestroyArray(iKeCoarser_dense); */
    }

    /* Free allocated memory for sK at the end */
    mxFree(sK);
}

/* AccumulationSmall function definition */
void AccumulationSmall(double *sK, double *tmpK, int *localMapping, mwSize nPerParentEleDOFs, mwSize perParentEleSons)
{
    /* Zero-initialize tmpK */
/*     for (int i = 0; i < nPerParentEleDOFs * nPerParentEleDOFs; i++) {
        tmpK[i] = 0.0;
    } */
	
    /* Perform the accumulation */
    for (int j = 0; j < perParentEleSons; j++) {
        for (mwSize k = 0; k < 576; k++) {
            int idx = localMapping[k + j * 576] - 1;	
            if (idx >= 0) {
				tmpK[idx] += sK[k + j * 576];
            }
        }
    }

    /* The result in tmpK is already reshaped implicitly, as MATLAB's reshape is just a reorganization of the memory */
}



/* Step 1: Compute B = interpolatingKe' * tmpK */
void ComputeIntermediateResult(double *B, const mwIndex *rowIndices, const mwIndex *colPointers, 
                               const double *values, const double *tmpK, mwSize nPerParentEleDOFs, mwSize nRows)
{
    /* Initialize B to zero manually */
/*     for (mwSize i = 0; i < nRows * nPerParentEleDOFs; i++) {
        B[i] = 0.0;
    } */
	
	for (mwSize col = 0; col < nRows; col++) {
        for (mwIndex k = colPointers[col]; k < colPointers[col + 1]; k++) {
            mwSize i = rowIndices[k];
            double val = values[k];
            for (mwSize row = 0; row < nPerParentEleDOFs; row++) {
                B[col * nPerParentEleDOFs + row] += val * tmpK[row * nPerParentEleDOFs + i];
            }
        }
    }
}

/* Step 2: Compute iKeCoarser = B * interpolatingKe */
void ComputeFinalResult(double *iKeCoarser, const double *B, const mwIndex *rowIndices, 
                        const mwIndex *colPointers, const double *values, mwSize nPerParentEleDOFs, mwSize nRows)
{
/*     for (mwSize i = 0; i < nRows * nRows; i++) {
        iKeCoarser[i] = 0.0;
    } */
    for (mwSize col = 0; col < nRows; col++) {
        for (mwIndex k = colPointers[col]; k < colPointers[col + 1]; k++) {
            mwSize i = rowIndices[k];
            double val = values[k];
            for (mwSize row = 0; row < nRows; row++) {
                iKeCoarser[col * nRows + row] += B[row * nPerParentEleDOFs + i] * val;
            }
        }
    }
}
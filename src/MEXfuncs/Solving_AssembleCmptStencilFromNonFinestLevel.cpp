#include "mex.h"
#include "matrix.h"
#include <omp.h>
#include <vector> 
#include <algorithm> 
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp /std:c++20" Solving_AssembleCmptStencilFromNonFinestLevel.cpp"
/* Main MEX function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Get input arguments */
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "5 inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "1 output required.");
    }

	// Retrieve the inputs
	const mxArray *KsPrevious_mx = prhs[0];
	const mxArray *elementUpwardMap_mx = prhs[1];
	const mxArray *interpolatingKe_mx = prhs[2]; 
	const mxArray *localMapping_mx = prhs[3];
	const mwSize numProjectNodes = static_cast<mwSize>(mxGetScalar(prhs[4]));
	
	// Basic constants
	const mwSize perEleDOFs = 24; // 24 DOFs per element (8 nodes * 3 DOFs)
	const mwSize nnzKe = perEleDOFs*perEleDOFs; //576 entries of Ke
	
    // Validate input dimensions
    if (!mxIsDouble(KsPrevious_mx)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:KsPreviousNotDouble",
                          "KsPrevious must be a double 24x24xN array.");
    }
    mwSize ndimsKs = mxGetNumberOfDimensions(KsPrevious_mx);
    if (ndimsKs != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:KsPreviousDim",
                          "KsPrevious must be 3D: 24x24xN.");
    }
    const mwSize *dimsKs = mxGetDimensions(KsPrevious_mx);
    if (dimsKs[0] != perEleDOFs || dimsKs[1] != perEleDOFs) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:KsPreviousSize",
                          "First two dimensions of KsPrevious must be 24x24.");
    }
    const mwSize nChildEles = dimsKs[2]; // number of fine-level elements/stencils	
    if (!mxIsInt32(elementUpwardMap_mx) || (mxGetN(elementUpwardMap_mx) != 8)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input elementUpwardMap must be an M*64 or M*8 int32 matrix.");
    }
    if (!mxIsSparse(interpolatingKe_mx) || !mxIsDouble(interpolatingKe_mx)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotSparse",  "Input interpolatingKe must be a sparse double matrix.");
    }	
    if (!mxIsInt32(localMapping_mx) || mxGetN(localMapping_mx) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input localMapping must be a column vector.");
    }
	
	// Get raw pointers
	double *KsPrevious = mxGetPr(KsPrevious_mx);
	int32_T *elementUpwardMap = (int32_T *)mxGetData(elementUpwardMap_mx);
	int32_T *localMapping = (int32_T *) mxGetData(localMapping_mx);	
	
    // Dimensions 
    mwSize nPerParentEleDOFs = 3 * numProjectNodes;
	const mwSize nnztmpK = nPerParentEleDOFs*nPerParentEleDOFs;
    mwSize perParentEleSons = mxGetN(elementUpwardMap_mx);  /* Number of columns in elementUpwardMap */
    mwSize nElesCurrent = mxGetM(elementUpwardMap_mx);  /* Number of rows in elementUpwardMap */

    /* Extract sparse matrix information */
    mwIndex *rowIndices = mxGetIr(interpolatingKe_mx);  /* Row indices of non-zero elements */
    mwIndex *colPointers = mxGetJc(interpolatingKe_mx); /* Column pointers for the sparse matrix */
    double *values = mxGetPr(interpolatingKe_mx);       /* Non-zero values in the sparse matrix */

    // Output: Ks (24 x 24 x nElesCurrent)
    mwSize dims[3] = {perEleDOFs, perEleDOFs, nElesCurrent};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *Ks = mxGetPr(plhs[0]);
    mwSize sK_size = nnzKe * perParentEleSons;
    
	#pragma omp parallel 
	{
		std::vector<double> sK(sK_size);
		std::vector<double> tmpK(nnztmpK);
		std::vector<double> B(perEleDOFs * nPerParentEleDOFs);
		std::vector<double> KeC(nnzKe);	
		double *sKptr  = sK.data();
		double *tmpKptr = tmpK.data();
		double *Bptr    = B.data();
		double *KeCptr  = KeC.data();
		
		#pragma omp for
		for (int i=0; i<(int)nElesCurrent; i++) {
			/* Zero-initialize sK before entering the inner loop */
			std::fill(sKptr, sKptr + sK_size, 0.0);
			/* Inner loop traverses the columns of elementUpwardMap */
			for (mwSize q=0; q<perParentEleSons; q++) {
				mwSize idx = elementUpwardMap[i + q*nElesCurrent]; 
				if (idx > 0) {
					for (mwSize col=0; col<perEleDOFs; col++) {
						for (mwSize row=0; row<perEleDOFs; row++)
						sKptr[col*perEleDOFs+row + q*nnzKe] = KsPrevious[col*perEleDOFs+row + (idx-1)*nnzKe];
					}
				}
			}
	
			/* Accumulate the values in sK using localMapping and reshape to sparse format */
			std::fill(tmpKptr, tmpKptr + nnztmpK, 0.0);
			for (mwSize j=0; j<perParentEleSons; j++) {
				for (mwSize k=0; k<nnzKe; k++) {
					mwSize idx = localMapping[k + j*nnzKe] - 1;
					tmpKptr[idx] += sKptr[k + j*nnzKe];
				}
			}
		
			// Compute transpose(interpolatingKe) * tmpK = B 
			std::fill(Bptr, Bptr + perEleDOFs * nPerParentEleDOFs, 0.0);
			for (mwSize col=0; col<perEleDOFs; col++) {
				for (mwIndex k=colPointers[col]; k<colPointers[col+1]; k++) {
					mwSize ss = rowIndices[k];
					double val = values[k];
					for (mwSize row=0; row<nPerParentEleDOFs; row++) {
						Bptr[col*nPerParentEleDOFs+row] += val * tmpKptr[row*nPerParentEleDOFs+ss];
					}
				}
			}
			// B * interpolatingKe = iKeCoarser
			std::fill(KeCptr, KeCptr+nnzKe, 0.0);			
			for (mwSize col=0; col<perEleDOFs; col++) {
				for (mwSize k=colPointers[col]; k<colPointers[col+1]; k++) {
					mwSize ss = rowIndices[k];
					double val = values[k];
					for (mwSize row=0; row<perEleDOFs; row++) {
						KeCptr[col*perEleDOFs+row] += Bptr[row*nPerParentEleDOFs+ss] * val;
					}
				}
			}
		
			/* Copy the dense result into the i-th page of Ks */
			for (mwSize row=0; row<perEleDOFs; row++) {
				for (mwSize col=0; col<perEleDOFs; col++) {
					Ks[row + col*perEleDOFs + i * nnzKe] = KeCptr[row + col * perEleDOFs];
				}
			}
		}
	}
}

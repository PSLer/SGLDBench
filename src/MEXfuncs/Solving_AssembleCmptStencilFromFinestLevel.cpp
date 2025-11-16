#include "mex.h"
#include "matrix.h"
#include <omp.h>
#include <vector>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_AssembleCmptStencilFromFinestLevel.cpp"
/*
	Compiling Command
	Windows: mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp" Solving_AssembleCmptStencilFromFinestLevel.cpp
	Linux: mex -v CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" Solving_AssembleCmptStencilFromFinestLevel.cpp
*/
/* Main MEX function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "6 inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "1 output required.");
    }
	
	// Retrieve the inputs
	const mxArray *iKe_mx = prhs[0];
	const mxArray *E_mx = prhs[1];
	const mxArray *elementUpwardMap_mx = prhs[2];
	const mxArray *interpolatingKe_mx = prhs[3]; 
	const mxArray *localMapping_mx = prhs[4];
	const mwSize numProjectNodes = static_cast<mwSize>(mxGetScalar(prhs[5]));
	
	// Basic constants
	const mwSize perEleDOFs = 24; // 24 DOFs per element (8 nodes * 3 DOFs)
	const mwSize nnzKe = perEleDOFs*perEleDOFs; //576 entries of Ke
	
    // Validate input dimensions
    if (!mxIsDouble(iKe_mx) || mxGetM(iKe_mx) != 24 || mxGetN(iKe_mx) != 24) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNot24x24", "Input Ke must be a 24x24 double matrix.");
    }
    if (!mxIsInt32(elementUpwardMap_mx) || (mxGetN(elementUpwardMap_mx) != 64 && mxGetN(elementUpwardMap_mx) != 8)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input elementUpwardMap must be an M*64 or M*8 int32 matrix.");
    }
    if (!mxIsSparse(interpolatingKe_mx) || !mxIsDouble(interpolatingKe_mx)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotSparse",  "Input interpolatingKe must be a sparse double matrix.");
    }	
    if (!mxIsInt32(localMapping_mx) || mxGetN(localMapping_mx) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input localMapping must be a column vector.");
    }
	
	// Get raw pointers
	double *iKe = mxGetPr(iKe_mx);
	double *eleModulus = mxGetPr(E_mx);
	int32_T *elementUpwardMap = (int32_T *)mxGetData(elementUpwardMap_mx);
	int32_T *localMapping = (int32_T *) mxGetData(localMapping_mx);	
	
    // Dimensions 
    mwSize nPerParentEleDOFs = 3 * numProjectNodes;
	const mwSize nnztmpK = nPerParentEleDOFs*nPerParentEleDOFs;
    mwSize perParentEleSons = mxGetN(elementUpwardMap_mx);  /* Number of columns in elementUpwardMap */
    mwSize nElesCurrent = mxGetM(elementUpwardMap_mx);  /* Number of rows in elementUpwardMap */	
	
    // Extract sparse matrix information for interpolatingKe
    mwIndex *rowIndices = mxGetIr(interpolatingKe_mx);  /* Row indices of non-zero elements */
    mwIndex *colPointers = mxGetJc(interpolatingKe_mx); /* Column pointers for the sparse matrix */
    double *values = mxGetPr(interpolatingKe_mx);       /* Non-zero values in the sparse matrix */

    // Output: Ks (24 x 24 x nElesCurrent)
    mwSize dims[3] = {perEleDOFs, perEleDOFs, nElesCurrent};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *Ks = mxGetPr(plhs[0]);
	const mwSize KsPageSize = nnzKe;
	
	// precompute per-child templates H_q (24x24 each) ----
    const mwSize H_blockSize = nnzKe; // 24*24
	std::vector<double> H(perParentEleSons * H_blockSize);
	// Temporary arrays for template construction (single-threaded)
    const mwSize sK_size = nnzKe * perParentEleSons; // 24*24*Q
	
    std::vector<double> sK_template(sK_size);
    std::vector<double> tmpK_template(nnztmpK);
    std::vector<double> B_template(perEleDOFs * nPerParentEleDOFs);
    std::vector<double> KeC_template(H_blockSize);	

    double *sK_t = sK_template.data();
    double *tmpK_t = tmpK_template.data();
    double *B_t = B_template.data();
    double *KeC_t = KeC_template.data();
	// Precompute templates H_q
	for (mwSize q=0; q<perParentEleSons; ++q) {
        std::fill(sK_t, sK_t + sK_size, 0.0);
        for (mwSize k=0; k<nnzKe; ++k) {
            sK_t[k + q*nnzKe] = iKe[k];
        }

        // tmpK is stored in column-major as nPerParentEleDOFs x nPerParentEleDOFs
        std::fill(tmpK_t, tmpK_t + nnztmpK, 0.0);
        for (mwSize j=0; j<perParentEleSons; ++j) {
            for (mwSize k=0; k<nnzKe; ++k) {
                int32_T lm = localMapping[k + j * nnzKe];
                if (lm<=0) continue; // safety
                mwSize flat = static_cast<mwSize>(lm - 1); // 0-based index
                tmpK_t[flat] += sK_t[k + j * nnzKe];
            }
        }

        // B(col,row) in column-major: B_t[col + row*24]
        std::fill(B_t, B_t + perEleDOFs * nPerParentEleDOFs, 0.0);
        for (mwSize col=0; col<perEleDOFs; ++col) {
            for (mwIndex kk=colPointers[col]; kk<colPointers[col+1]; ++kk) {
                mwSize r = rowIndices[kk]; // row index in [0..nPerParentEleDOFs-1]
                double val = values[kk];   // R(r, col)
                // B(col, row) += R(r,col) * tmpK(row, r)
                for (mwSize row=0; row<nPerParentEleDOFs; ++row) {
                    // tmpK(row, r) = tmpK_t[row + r * nPerParentEleDOFs]
                    B_t[col + row*perEleDOFs] += val * tmpK_t[row + r * nPerParentEleDOFs];
                }
            }
        }
		
        // KeC(row,col) in column-major: KeC_t[row + col*24]
        std::fill(KeC_t, KeC_t + H_blockSize, 0.0);
        for (mwSize colR=0; colR<perEleDOFs; ++colR) {
            for (mwIndex kk=colPointers[colR]; kk<colPointers[colR+1]; ++kk) {
                mwSize r = rowIndices[kk];
                double val = values[kk];   
                for (mwSize row=0; row<perEleDOFs; ++row) {
                    KeC_t[row + colR*perEleDOFs] += B_t[row + r*perEleDOFs] * val;
                }
            }
        }
		
        double *Hq = &H[q * H_blockSize];
        std::copy(KeC_t, KeC_t + H_blockSize, Hq);		
	}
	
    // Main parallel loop: linear combination of templates
	#pragma omp parallel 
	{
		std::vector<double> KeC(H_blockSize);
        double *KeCptr = KeC.data();
		
		 #pragma omp for schedule(static)
		for (int i=0; i < (int)nElesCurrent; i++) {
			// Zero parent coarse stiffness
            std::fill(KeCptr, KeCptr + H_blockSize, 0.0);
	
			// Accumulate contributions from each child position q
            for (mwSize q=0; q<perParentEleSons; ++q) {
                int32_T idxChild = elementUpwardMap[i + q*nElesCurrent];  // MATLAB column-major

                if (idxChild==0) {
                    continue; // ghost or void child
                }

                // idxChild is 1-based index into eleModulus
                double E_value = eleModulus[idxChild-1];
                const double *Hq = &H[q * H_blockSize];

                // KeC += E_value * Hq
                for (mwSize t=0; t<H_blockSize; ++t) {
                    KeCptr[t] += E_value * Hq[t];
                }
            }
			
            // Copy into Ks(:,:,i)
            // Ks is 24 x 24 x nElesCurrent in column-major:
            // Ks(row, col, i) = Ks[row + col*24 + i*24*24]
            double *KsPage = &Ks[(mwSize)i * KsPageSize];
            std::copy(KeCptr, KeCptr+H_blockSize, KsPage);
		}
	}	
}

#include "mex.h"
#include <omp.h>
#include <stdlib.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp:experimental /std:c++20" Solving_KbyU_MatrixFree8x8_mex.cpp"
inline void apply_Ke_Ue_8(const double* __restrict Ke, const double* __restrict Ue, double* __restrict Ye);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "Five inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "One output required.");
    }

    // Retrieve the inputs
    const mxArray *U_mx        = prhs[0];
    const mxArray *eNodMat_mx  = prhs[1];
    const mxArray *Ke_mx       = prhs[2];
    const mxArray *colorArray  = prhs[3];
	
    double *U = mxGetPr(U_mx);
    int32_T *eNodMat = (int32_T *)mxGetData(eNodMat_mx);
    double *Ke = mxGetPr(Ke_mx);

    mwSize numDOFs = mxGetM(U_mx);
    mwSize numElements = mxGetM(eNodMat_mx);
	mwSize nColors = mxGetNumberOfElements(colorArray);

    // Validate input dimensions
    if (!mxIsDouble(U_mx) || mxGetN(U_mx) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input U must be a column vector.");
    }
    if (!mxIsInt32(eNodMat_mx) || mxGetN(eNodMat_mx) != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input eNodMat must be an Mx8 int32 matrix.");
    }
    if (!mxIsDouble(Ke_mx) || mxGetM(Ke_mx) != 8 || mxGetN(Ke_mx) != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNot24x24", "Input Ke must be a 24x24 double matrix.");
    }
    if (!mxIsCell(colorArray)) {
        mexErrMsgIdAndTxt("Solving_KbyU:colorArray", "colorArray must be a cell array.");
    }	
	

    // Create output array Y with dimensions [numDOFs, 1]
    plhs[0] = mxCreateDoubleMatrix(numDOFs, 1, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    // Initialize Y to zeros
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < numDOFs; ++i) {
        Y[i] = 0.0;
    }
	
	for (mwSize icc=0; icc<nColors; ++icc) {
		const mxArray *elesUnique = mxGetCell(colorArray, icc);
        if (elesUnique == nullptr) continue;
        if (!mxIsInt32(elesUnique) || mxIsComplex(elesUnique)) {
            mexErrMsgIdAndTxt("Solving_KbyU:colorCell", "Each color cell must be a real int32 vector.");
        }
		
		int32_T *elemIdx = (int32_T *)mxGetData(elesUnique);
		mwSize iNumElesUnique = mxGetNumberOfElements(elesUnique);	
		
		#pragma omp parallel for schedule(static)
		for (int ieLocal=0; ieLocal<iNumElesUnique; ++ieLocal) {
            // Local arrays per element (on stack, thread-private)
            int32_T ieGlobal = elemIdx[ieLocal]-1;
			int32_T ieNod[8];
            double  Ue[8];
            double  Ye[8];
			for (int j = 0; j < 8; ++j) {
				ieNod[j] = eNodMat[ieGlobal + j * numElements];
			}
				
			// Create Ue matrix with dimensions [24,1]
			for (int j = 0; j < 8; ++j) {
                mwIndex node = (mwIndex)ieNod[j] - 1;
				Ue[j] = U[node];
			}
			
			// Create Ye matrix with dimensions [24,1] Ye = Ke*Ue
			apply_Ke_Ue_8(Ke, Ue, Ye);

			// Accumulate the elements in Ye into Y
            for (int j = 0; j < 8; ++j) {
                mwIndex node = (mwIndex)ieNod[j] - 1;
				Y[node] += Ye[j];
            }			
		}
	}
}

inline void apply_Ke_Ue_8(const double* __restrict Ke, const double* __restrict Ue, double* __restrict Ye)
{
	if (0) {
		for (int j = 0; j < 8; ++j) {
			double sum = 0.0;
			for (int k = 0; k < 8; ++k) {
				sum += Ue[k] * Ke[k * 8 + j];
			}
			Ye[j] = sum;
		}		
	} 
	else {
		// initialize Ye
		for (int j = 0; j < 8; ++j) {
			Ye[j] = 0.0;
		}
	
		// Ke * Ue
		for (int k = 0; k < 8; ++k) {
			double uk = Ue[k];
			const double* Kcol = &Ke[k * 8];
	
			#pragma omp simd
			for (int j = 0; j < 8; ++j) {
				Ye[j] += Kcol[j] * uk;
			}
		}
	}
}

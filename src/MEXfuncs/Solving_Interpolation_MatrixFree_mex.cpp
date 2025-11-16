#include "mex.h"
#include <stdint.h>
#include <omp.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp:experimental /std:c++20" Solving_Interpolation_MatrixFree_mex.cpp"
inline void apply_perEleInterpolation(const double* __restrict RI, const int numChildNodesPerEle, const double* __restrict UeCoarse, double* __restrict UeFiner);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 9) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "9 inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "1 output required.");
    }
    // Retrieve the inputs
    const mxArray *xCoarser_mx = prhs[0];
	const mxArray *transferMat_mx = prhs[1];
	const mxArray *RI_mx = prhs[2];
	const mxArray *eNodMat_mx = prhs[3];
	const mxArray *solidNodeMapCoarser2Finer_mx = prhs[6];
	const mxArray *transferMatCoeffi_mx = prhs[7];
	const mxArray *colorArray  = prhs[8];
	
	double *xCoarser = mxGetPr(xCoarser_mx);
	int32_T *transferMat = (int32_T *)mxGetData(transferMat_mx);
	double *RI = mxGetPr(RI_mx);
	int32_T *eNodMat = (int32_T *)mxGetData(eNodMat_mx);
	mwSize numNodesFine = (mwSize)mxGetScalar(prhs[4]);
	mwSize numNodesFineGhost = (mwSize)mxGetScalar(prhs[5]);
	int32_T *solidNodeMapCoarser2Finer = (int32_T *)mxGetData(solidNodeMapCoarser2Finer_mx);
	double *transferMatCoeffi = mxGetPr(transferMatCoeffi_mx);
	
    // Validate input dimensions
    if (!mxIsDouble(xCoarser_mx) || mxGetN(xCoarser_mx) != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input xCoarser must be a n-by-3 matrix.");
    }
    if (!mxIsInt32(transferMat_mx) || (mxGetM(transferMat_mx) != 125 && mxGetM(transferMat_mx) != 27)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input transferMat must be an 27xM or 125xM int32 matrix.");
    }
    if (!mxIsDouble(RI_mx) || mxGetN(RI_mx) != 8 || (mxGetM(RI_mx) != 125 && mxGetM(RI_mx) != 27)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMx8", "Input RI must be a 125/27 x8 double matrix.");
    }
    if (!mxIsInt32(eNodMat_mx) || mxGetN(eNodMat_mx) != 8) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input eNodMat must be an Mx8 int32 matrix.");
    }
    if (!mxIsInt32(solidNodeMapCoarser2Finer_mx) || mxGetN(solidNodeMapCoarser2Finer_mx) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input solidNodeMapCoarser2Finer must be a column vector.");
    }	
    if (!mxIsDouble(transferMatCoeffi_mx) || mxGetN(transferMatCoeffi_mx) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input transferMatCoeffi must be a column vector.");
    }
    if (!mxIsCell(colorArray)) {
        mexErrMsgIdAndTxt("Solving_KbyU:colorArray", "colorArray must be a cell array.");
    }		
	
	mwSize numDOFsPerNode = mxGetN(xCoarser_mx); //constant 3
	mwSize numNodesCoarse = mxGetM(xCoarser_mx);
	mwSize numElesCoarse = mxGetM(eNodMat_mx); // Number of rowsTransferMat in transferMat
	mwSize nColors = mxGetNumberOfElements(colorArray);
	mwSize mm = mxGetM(RI_mx); // Number of rowsTransferMat in multiGridOperatorRI (125 or 27)
	int numChildNodesPerEle = static_cast<int>(mm);
	
	//default output
	plhs[0] = mxCreateDoubleMatrix(numNodesFine, numDOFsPerNode, mxREAL);
	double *Y = mxGetPr(plhs[0]);

#if 1 //efficiency-prior
	double *YGhost = (double *)mxMalloc(numNodesFineGhost*numDOFsPerNode * sizeof(double));
	int numDOFsFineGhost = static_cast<int>(numNodesFineGhost*numDOFsPerNode);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < numDOFsFineGhost; ++i) {
        YGhost[i] = 0.0;
    }
	
	// Interpolating
	for (mwSize icc=0; icc<nColors; ++icc) { //per Color
		const mxArray *elesUnique = mxGetCell(colorArray, icc);
		if (elesUnique == nullptr) continue;
        if (!mxIsInt32(elesUnique) || mxIsComplex(elesUnique)) {
            mexErrMsgIdAndTxt("Solving_Interpolation:colorCell", "Each color cell must be a real int32 vector.");
        }
		int32_T *elemIdx = (int32_T *)mxGetData(elesUnique);
		mwSize iNumElesUnique = mxGetNumberOfElements(elesUnique);
		
		//Initialize child variables
		#pragma omp parallel for schedule(static)
		for (int ieLocal=0; ieLocal<iNumElesUnique; ++ieLocal) {
            int32_T ieGlobal = elemIdx[ieLocal]-1;
			int32_T ieNodCoarse[8];
			double UeX[8], UeY[8], UeZ[8];
			for (int j=0; j<8; ++j) {
				ieNodCoarse[j] = eNodMat[ieGlobal + j*numElesCoarse];
			}
			
			// Gather element-wise vector at coarser level
			for (int j=0; j<8; ++j) {
				mwIndex jNode = (mwIndex)ieNodCoarse[j] - 1;
				UeX[j] = xCoarser[jNode];
				UeY[j] = xCoarser[jNode+numNodesCoarse];
				UeZ[j] = xCoarser[jNode+2*numNodesCoarse];					
			}
			
			// Conduct element-wise interpolation
			double YeX[125], YeY[125], YeZ[125];
			apply_perEleInterpolation(RI, numChildNodesPerEle, UeX, YeX);
			apply_perEleInterpolation(RI, numChildNodesPerEle, UeY, YeY);
			apply_perEleInterpolation(RI, numChildNodesPerEle, UeZ, YeZ);
			
			
			int32_T ieNodFine[125];
			for (int j=0; j<numChildNodesPerEle; ++j) {
				ieNodFine[j] = transferMat[numChildNodesPerEle*ieGlobal+j];
			}			
			for (int j=0; j<numChildNodesPerEle; ++j) {
				mwIndex jNode = (mwIndex)ieNodFine[j] - 1;
				YGhost[jNode] += YeX[j];
				YGhost[jNode+numNodesFineGhost] += YeY[j];
				YGhost[jNode+2*numNodesFineGhost] += YeZ[j];				
			}
		}		
	}
	for (int ss=0; ss<numDOFsPerNode; ++ss) {
		mwIndex iDOFghost = (mwIndex)ss*numNodesFineGhost;
		mwIndex iDOF = (mwIndex)ss*numNodesFine;
		#pragma omp parallel for
		for (int j=0; j<numNodesFine; ++j) {
			mwIndex jDOFghost = solidNodeMapCoarser2Finer[j]-1;
			Y[j+iDOF] = YGhost[jDOFghost + iDOFghost] / transferMatCoeffi[jDOFghost];
		}
	}
#else //memory-prior
	double *YGhost = (double *)mxMalloc(numNodesFineGhost * sizeof(double));
	int num1DOFFineGhost = static_cast<int>(numNodesFineGhost);
	
	// Interpolating
	for (int ss=0; ss<numDOFsPerNode; ++ss) {
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < num1DOFFineGhost; ++i) {
			YGhost[i] = 0.0;
		}		
		int sDOFcoarse = ss*numNodesCoarse;
		int sDOFfineGhost = ss*numNodesFineGhost;
		for (mwSize icc=0; icc<nColors; ++icc) {
			const mxArray *elesUnique = mxGetCell(colorArray, icc);
			if (elesUnique == nullptr) continue;
			if (!mxIsInt32(elesUnique) || mxIsComplex(elesUnique)) {
				mexErrMsgIdAndTxt("Solving_Interpolation:colorCell", "Each color cell must be a real int32 vector.");
			}
			int32_T *elemIdx = (int32_T *)mxGetData(elesUnique);
			mwSize iNumElesUnique = mxGetNumberOfElements(elesUnique);
			
			//Initialize child variables
			#pragma omp parallel for schedule(static)
			for (int ieLocal=0; ieLocal<iNumElesUnique; ++ieLocal) {
				int32_T ieGlobal = elemIdx[ieLocal]-1;
				int32_T ieNodCoarse[8];
				double Ue[8];
				for (int j=0; j<8; ++j) {
					ieNodCoarse[j] = eNodMat[ieGlobal + j*numElesCoarse];
				}
				
				// Gather element-wise vector at coarser level
				for (int j=0; j<8; ++j) {
					mwIndex jNode = (mwIndex)ieNodCoarse[j] - 1;
					Ue[j] = xCoarser[jNode+sDOFcoarse];				
				}
				
				// Conduct element-wise interpolation
				double Ye[125];
				apply_perEleInterpolation(RI, numChildNodesPerEle, Ue, Ye);		
				
				int32_T ieNodFine[125];
				for (int j=0; j<numChildNodesPerEle; ++j) {
					ieNodFine[j] = transferMat[numChildNodesPerEle*ieGlobal+j];
				}			
				for (int j=0; j<numChildNodesPerEle; ++j) {
					mwIndex jNode = (mwIndex)ieNodFine[j] - 1;
					YGhost[jNode] += Ye[j];			
				}
			}			
		}
		
		mwIndex iDOF = (mwIndex)ss*numNodesFine;
		#pragma omp parallel for
		for (int j = 0; j < numNodesFine; ++j) {
			mwIndex jDOFghost = solidNodeMapCoarser2Finer[j]-1;
			Y[j+iDOF] = YGhost[jDOFghost] / transferMatCoeffi[jDOFghost];
		}		
	}
#endif	
	mxFree(YGhost);
}

inline void apply_perEleInterpolation(const double* __restrict RI, const int numChildNodesPerEle, const double* __restrict UeCoarse, double* __restrict UeFiner) {
	if (0) {
		for (int j=0; j<numChildNodesPerEle; ++j) {
			double sum = 0.0;
			for (int k=0; k<8; ++k) {
				sum += UeCoarse[k] * RI[k * numChildNodesPerEle + j];
			}
			UeFiner[j] = sum;
		}		
	} 
	else {
		// initialize UeFiner
		for (int j=0; j<numChildNodesPerEle; ++j) {
			UeFiner[j] = 0.0;
		}
	
		// RI * UeCoarse
		for (int k=0; k<8; ++k) {
			double uk = UeCoarse[k];
			const double* RIcol = &RI[k * numChildNodesPerEle];
	
			#pragma omp simd
			for (int j=0; j<numChildNodesPerEle; ++j) {
				UeFiner[j] += RIcol[j] * uk;
			}
		}
	}	
}
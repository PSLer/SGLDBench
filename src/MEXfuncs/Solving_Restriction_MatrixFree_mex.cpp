#include "mex.h"
#include <stdint.h>
#include <omp.h>
// Compile with "mex -largeArrayDims COMPFLAGS="$COMPFLAGS /openmp:experimental /std:c++20" Solving_Restriction_MatrixFree_mex.cpp"
inline void apply_perEleRestriction(const double* __restrict RIT, const int numChildNodesPerEle, const double* __restrict UeFiner, double* __restrict UeCoarse);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for the correct number of input arguments
    if (nrhs != 9) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nrhs", "9 inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:nlhs", "1 output required.");
    }
    // Retrieve the inputs
    const mxArray *rFiner_mx = prhs[0];
	const mxArray *transferMat_mx = prhs[1];
	const mxArray *RIT_mx = prhs[2];
	const mxArray *eNodMat_mx = prhs[3];
	const mxArray *solidNodeMapCoarser2Finer_mx = prhs[6];
	const mxArray *transferMatCoeffi_mx = prhs[7];
	const mxArray *colorArray  = prhs[8];

	double *rFiner = mxGetPr(rFiner_mx);
	int32_T *transferMat = (int32_T *)mxGetData(transferMat_mx);
	double *RIT = mxGetPr(RIT_mx);
	int32_T *eNodMat = (int32_T *)mxGetData(eNodMat_mx);
	mwSize numNodesCoarse = (mwSize)mxGetScalar(prhs[4]);
	mwSize numNodesFineGhost = (mwSize)mxGetScalar(prhs[5]);
	int32_T *solidNodeMapCoarser2Finer = (int32_T *)mxGetData(solidNodeMapCoarser2Finer_mx);
	double *transferMatCoeffi = mxGetPr(transferMatCoeffi_mx);

    // Validate input dimensions
    if (!mxIsDouble(rFiner_mx) || mxGetN(rFiner_mx) != 3) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotVector", "Input rFiner must be a n-by-3 matrix.");
    }
    if (!mxIsInt32(transferMat_mx) || (mxGetM(transferMat_mx) != 125 && mxGetM(transferMat_mx) != 27)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMatrix", "Input transferMat must be an 27xM or 125xM int32 matrix.");
    }
    if (!mxIsDouble(RIT_mx) || mxGetM(RIT_mx) != 8 || (mxGetN(RIT_mx) != 125 && mxGetN(RIT_mx) != 27)) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotMx8", "Input RIT must be a 8x 125/27 double matrix.");
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

	mwSize numDOFsPerNode = mxGetN(rFiner_mx); //constant 3
	mwSize numNodesFine = mxGetM(rFiner_mx);
	mwSize numElesCoarse = mxGetM(eNodMat_mx); // Number of rowsTransferMat in transferMat
	mwSize nColors = mxGetNumberOfElements(colorArray);
	mwSize mm = mxGetN(RIT_mx); // (125 or 27)
	int numChildNodesPerEle = static_cast<int>(mm);
	int numDOFsCoarse = static_cast<int>(numNodesCoarse*numDOFsPerNode);
	
	//default output
	plhs[0] = mxCreateDoubleMatrix(numNodesCoarse, numDOFsPerNode, mxREAL);
	double *Y = mxGetPr(plhs[0]);
	#pragma omp parallel for
	for (int i=0; i<numDOFsCoarse; ++i) {
		Y[i]=0;
	}
	
#if 1 //efficiency-prior
	double *rFinerGhost = (double *)mxMalloc(numNodesFineGhost*numDOFsPerNode * sizeof(double));
	int numDOFsFineGhost = static_cast<int>(numNodesFineGhost*numDOFsPerNode);
	#pragma omp parallel for
	for (int i=0; i<numDOFsFineGhost; ++i) {
		rFinerGhost[i] = 0.0;
	}
	for (int ss=0; ss<numDOFsPerNode; ++ss) {
		mwIndex iDOFghost = (mwIndex)ss*numNodesFineGhost;
		mwIndex iDOF = (mwIndex)ss*numNodesFine;
		#pragma omp parallel for
		for (int j=0; j<numNodesFine; ++j) {
			mwIndex jNode = solidNodeMapCoarser2Finer[j]-1;
			rFinerGhost[jNode+iDOFghost] = rFiner[iDOF+j]/transferMatCoeffi[jNode];
		}		
	}
	
	// Restriction
	for (mwSize icc=0; icc<nColors; ++icc) { //per Color
		const mxArray *elesUnique = mxGetCell(colorArray, icc);
		if (elesUnique == nullptr) continue;
        if (!mxIsInt32(elesUnique) || mxIsComplex(elesUnique)) {
            mexErrMsgIdAndTxt("Solving_Interpolation:colorCell", "Each color cell must be a real int32 vector.");
        }
		int32_T *elemIdx = (int32_T *)mxGetData(elesUnique);
		mwSize iNumElesUnique = mxGetNumberOfElements(elesUnique);
		
		#pragma omp parallel for schedule(static)
		for (int ieLocal=0; ieLocal<iNumElesUnique; ++ieLocal) {
			int32_T ieGlobal = elemIdx[ieLocal]-1;
			int32_T ieNodFine[125];
			double UeX[125], UeY[125], UeZ[125];
			for (int j=0; j<numChildNodesPerEle; ++j) {
				ieNodFine[j] = transferMat[numChildNodesPerEle*ieGlobal+j];
			}
			
			// Gather element-wise vector at finer level
			for (int j=0; j<numChildNodesPerEle; ++j) {
				mwIndex jNode = (mwIndex)ieNodFine[j] - 1;
				UeX[j] = rFinerGhost[jNode];
				UeY[j] = rFinerGhost[jNode+numNodesFineGhost];
				UeZ[j] = rFinerGhost[jNode+2*numNodesFineGhost];				
			}

			// Conduct element-wise restriction
			double YeX[8], YeY[8], YeZ[8];
			apply_perEleRestriction(RIT, numChildNodesPerEle, UeX, YeX);
			apply_perEleRestriction(RIT, numChildNodesPerEle, UeY, YeY);
			apply_perEleRestriction(RIT, numChildNodesPerEle, UeZ, YeZ);
			
			int32_T ieNodCoarse[8];
			for (int j=0; j<8; ++j) {
				ieNodCoarse[j] = eNodMat[ieGlobal + j*numElesCoarse];
			}
			for (int j=0; j<8; ++j) {
				mwIndex jNode = (mwIndex)ieNodCoarse[j] - 1;
				Y[jNode] += YeX[j];
				Y[jNode+numNodesCoarse] += YeY[j];
				Y[jNode+2*numNodesCoarse] += YeZ[j];				
			}			
		}
	}
#else //memory-prior
	double *rFinerGhost = (double *)mxMalloc(numNodesFineGhost * sizeof(double));
	int num1DOFFineGhost = static_cast<int>(numNodesFineGhost);
	for (int ss=0; ss<numDOFsPerNode; ++ss) {
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < num1DOFFineGhost; ++i) {
			rFinerGhost[i] = 0.0;
		}
		mwIndex ssDOFcoarse = (mwIndex)ss*numNodesCoarse;
		mwIndex ssDOFiner = (mwIndex)ss*numNodesFine;
		
		#pragma omp parallel for
		for (int j=0; j<numNodesFine; ++j) {
			mwIndex jNode = solidNodeMapCoarser2Finer[j]-1;
			rFinerGhost[jNode] = rFiner[ssDOFiner+j]/transferMatCoeffi[jNode];
		}
		// Restriction
		for (mwSize icc=0; icc<nColors; ++icc) { //per Color
			const mxArray *elesUnique = mxGetCell(colorArray, icc);
			if (elesUnique == nullptr) continue;
			if (!mxIsInt32(elesUnique) || mxIsComplex(elesUnique)) {
				mexErrMsgIdAndTxt("Solving_Interpolation:colorCell", "Each color cell must be a real int32 vector.");
			}
			int32_T *elemIdx = (int32_T *)mxGetData(elesUnique);
			mwSize iNumElesUnique = mxGetNumberOfElements(elesUnique);

			#pragma omp parallel for schedule(static)
			for (int ieLocal=0; ieLocal<iNumElesUnique; ++ieLocal) {
				int32_T ieGlobal = elemIdx[ieLocal]-1;
				int32_T ieNodFine[125];
				double Ue[125];
				for (int j=0; j<numChildNodesPerEle; ++j) {
					ieNodFine[j] = transferMat[numChildNodesPerEle*ieGlobal+j];
				}
				
				// Gather element-wise vector at finer level
				for (int j=0; j<numChildNodesPerEle; ++j) {
					mwIndex jNode = (mwIndex)ieNodFine[j] - 1;
					Ue[j] = rFinerGhost[jNode];				
				}
	
				// Conduct element-wise restriction
				double Ye[8];
				apply_perEleRestriction(RIT, numChildNodesPerEle, Ue, Ye);
				
				int32_T ieNodCoarse[8];
				for (int j=0; j<8; ++j) {
					ieNodCoarse[j] = eNodMat[ieGlobal + j*numElesCoarse];
				}
				for (int j=0; j<8; ++j) {
					mwIndex jNode = (mwIndex)ieNodCoarse[j] - 1;
					Y[jNode+ssDOFcoarse] += Ye[j];				
				}			
			}
		}
	}
#endif
	mxFree(rFinerGhost);
}

inline void apply_perEleRestriction(const double* __restrict RIT, const int numChildNodesPerEle, const double* __restrict UeFiner, double* __restrict UeCoarse) {
	if (0) {
		for (int j=0; j<8; ++j) {
			double sum = 0.0;
			for (int k=0; k<numChildNodesPerEle; ++k) {
				sum += UeFiner[k] * RIT[k * 8 + j];
			}
			UeCoarse[j] = sum;
		}		
	} 
	else {
		// initialize UeCoarse
		for (int j=0; j<8; ++j) {
			UeCoarse[j] = 0.0;
		}
	
		// RIT * UeFiner
		for (int k=0; k<numChildNodesPerEle; ++k) {
			double uk = UeFiner[k];
			const double* RIcol = &RIT[k * 8];
	
			#pragma omp simd
			for (int j=0; j<8; ++j) {
				UeCoarse[j] += RIcol[j] * uk;
			}
		}
	}	
}
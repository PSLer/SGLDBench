#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Check the number of inputs and outputs
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:nrhs", "Five input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:nlhs", "Too many output arguments.");
    }

    // Get the input 3D matrix KsPrevious
    const mxArray *KsPrevious = prhs[0];
    if (!mxIsDouble(KsPrevious) || mxGetNumberOfDimensions(KsPrevious) != 3) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:inputNot3D", "KsPrevious must be a 3D double matrix.");
    }
    
    // Get dimensions of KsPrevious
    mwSize dims[3];
    mxGetDimensions(KsPrevious, dims);
    mwSize N = dims[2];
    
    // Get the number of elements and project nodes
    const mxArray *elementUpwardMap = prhs[1];
    const mxArray *interpolatingKe = prhs[2];
    const mxArray *localMapping = prhs[3];
    const mxArray *numProjectNodesScalar = prhs[4];
    
    if (!mxIsDouble(numProjectNodesScalar) || mxGetNumberOfElements(numProjectNodesScalar) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:numProjectNodes", "numProjectNodes must be a scalar.");
    }
    
    // Extract the number of project nodes
    int numProjectNodes = (int)mxGetScalar(numProjectNodesScalar);
    int nPerParentEleDOFs = 3 * numProjectNodes;
    
    // Check elementUpwardMap dimensions
    if (!mxIsInt32(elementUpwardMap) || mxGetNumberOfDimensions(elementUpwardMap) != 2) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:elementUpwardMap", "elementUpwardMap must be a 2D int32 matrix.");
    }
    mwSize* elementUpwardMapDims = mxGetDimensions(elementUpwardMap);
    mwSize nElesCurrent = elementUpwardMapDims[0];
    mwSize perParentEleSons = elementUpwardMapDims[1];
    
    // Check interpolatingKe dimensions
    if (!mxIsSparse(interpolatingKe) || mxGetM(interpolatingKe) != nPerParentEleDOFs || mxGetN(interpolatingKe) != 24) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:interpolatingKe", "interpolatingKe must be a sparse matrix of size nPerParentEleDOFs x 24.");
    }
    
    // Check localMapping dimensions
    if (!mxIsInt32(localMapping) || mxGetNumberOfDimensions(localMapping) != 2) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:localMapping", "localMapping must be a 2D int32 matrix.");
    }
    mwSize* localMappingDims = mxGetDimensions(localMapping);
    if (localMappingDims[0] != 24 * 24 * perParentEleSons || localMappingDims[1] != 1) {
        mexErrMsgIdAndTxt("MyToolbox:initializeKs:localMapping", "localMapping must be a column vector with dimensions (24*24*perParentEleSons) x 1.");
    }

    // Initialize the output 3D matrix Ks with dimensions 24x24xN
    mxArray *Ks = mxCreateNumericMatrix(24, 24, mxDOUBLE_CLASS, mxREAL);
    mwSize KsDims[3] = {24, 24, nElesCurrent};
    mxArray *KsArray = mxCreateNumericArray(3, KsDims, mxDOUBLE_CLASS, mxREAL);
    double *KsData = mxGetPr(KsArray);

    // Initialize Ks to zero
    mxSetPr(KsArray, KsData);
    for (mwSize i = 0; i < 24 * 24 * nElesCurrent; ++i) {
        KsData[i] = 0.0;
    }

    // Set the output
    plhs[0] = KsArray;
}

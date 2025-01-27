function vonMisesStressPerElement = FEA_ComputePerElementVonMisesStress(cartesianStressField)
	global meshHierarchy_;
	numElements = meshHierarchy_(1).numElements;
	
	domiDirField = [];
	% [cartesianStressField,~] = FEA_StressAnalysis();
	domiDirField = zeros(numElements,3);
	
	eNodMat = meshHierarchy_(1).eNodMat;
	shapeFuncsAtCentroid = FEA_ShapeFunction(0.0, 0.0, 0.0);
	cartesianStressFieldPerEle = zeros(numElements,6);
	
	if isempty(gcp('nocreate')), parpool('threads'); end		
	parfor ii=1:numElements
		iCartesianStressEleNodes = cartesianStressField(eNodMat(ii,:), :);
		cartesianStressFieldPerEle(ii,:) = shapeFuncsAtCentroid * iCartesianStressEleNodes;
	end
	vonMisesStressPerElement = FEA_ComputeVonMisesStress(cartesianStressFieldPerEle);
	vonMisesStressPerElement = vonMisesStressPerElement.^(1/2); %%Exaggerated display
end
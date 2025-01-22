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
	% vonMisesStressPerElement = sqrt(0.5*((cartesianStressFieldPerEle(:,1)-cartesianStressFieldPerEle(:,2)).^2 + ...
		% (cartesianStressFieldPerEle(:,2)-cartesianStressFieldPerEle(:,3)).^2 + (cartesianStressFieldPerEle(:,3)...
			% -cartesianStressFieldPerEle(:,1)).^2 ) + 3*( cartesianStressFieldPerEle(:,6).^2 + cartesianStressFieldPerEle(:,4).^2 + ...
				% cartesianStressFieldPerEle(:,5).^2 ));	
	vonMisesStressPerElement = FEA_ComputeVonMisesStress(cartesianStressFieldPerEle);
end
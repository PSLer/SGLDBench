function domiDirField = Common_ExtractDominantDirectionsFromPrincipalStressDirections(cartesianStressField)
	global meshHierarchy_;
	numElements = meshHierarchy_(1).numElements;
	
	domiDirField = [];
	% [cartesianStressField,~] = FEA_StressAnalysis();
	domiDirField = zeros(numElements,3);
	
	eNodMat = meshHierarchy_(1).eNodMat;
	shapeFuncsAtCentroid = FEA_ShapeFunction(0.0, 0.0, 0.0);
	principalStressFieldPerEle = zeros(numElements,12);
	
	if isempty(gcp('nocreate')), parpool('threads'); end		
	parfor ii=1:numElements
		iCartesianStressEleNodes = cartesianStressField(eNodMat(ii,:), :);
		iCartesianStressEle = shapeFuncsAtCentroid * iCartesianStressEleNodes;
		principalStressFieldPerEle(ii,:) = FEA_ComputePrincipalStress(iCartesianStressEle);
	end
	clear cartesianStressField	
	parfor ii=1:numElements
		iPS = principalStressFieldPerEle(ii,:);
		iPSamp = abs(iPS([1 5 9]));
		[~, whichPSisDominant] = max(iPSamp);
		switch whichPSisDominant
			case 1
				domiDirField(ii,:) = iPS([2 3 4]);
			case 2
				domiDirField(ii,:) = iPS([6 7 8]);
				disp('I am here!');
			case 3
				domiDirField(ii,:) = iPS([10 11 12]);
		end
	end
	clear principalStressFieldPerEle	
end
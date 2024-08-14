function domiDirField = Common_ExtractDominantDirectionsFromPrincipalStressDirections()
	global meshHierarchy_;
	global cartesianStressField_;
	numElements = meshHierarchy_(1).numElements;
	domiDirField = zeros(numElements,3);
	
	FEA_StressAnalysis();
	cartesianStressField = cartesianStressField_; clear -global cartesianStressField_
	eNodMat = meshHierarchy_(1).eNodMat;
	shapeFuncsAtCentroid = FEA_ShapeFunction(0.0, 0.0, 0.0);
	principalStressFieldPerEle = zeros(numElements,12);
	
	delete(gcp('nocreate')); p = parpool('threads', feature('numcores'));
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
				domiDirField(:,ii) = iPS([2 3 4]);
			case 2
				domiDirField(:,ii) = iPS([6 7 8]);
			case 3
				domiDirField(:,ii) = iPS([10 11 12]);
		end
	end
	delete(p);
	clear principalStressFieldPerEle	
end
function alignmentMetricVolume = Common_ComputeStressAlignmentDeviation(dominantDirSolid, dominantDirDesign)
	global meshHierarchy_;
	global densityLayout_;
	
	numElements = meshHierarchy_(1).numElements;
	alignmentMetric = zeros(size(dominantDirSolid,1),1);
	if size(dominantDirSolid,1)~=numElements || numElements~=size(dominantDirDesign,1)
		error('Un-matched Datasets!');
	end
	
	voidElementsInDesign = densityLayout_<0.1;
	delete(gcp('nocreate')); p = parpool('threads', feature('numcores'));
	parfor ii=1:numElements
		if ~voidElementsInDesign
			v1 = dominantDirSolid(ii,:);
			v2 = dominantDirDesign(ii,:);
			iDirCos = v1 * v2' / norm(v1) / norm(v2);
			incAngSin = sqrt(1 - iDirCos^2);
			alignmentMetric(ii,1) = incAngSin;
		end
	end
	delete(p);	

	alignmentMetricVolume = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	alignmentMetricVolume(meshHierarchy_(1).eleMapBack,1) = alignmentMetric;
	alignmentMetricVolume = reshape(alignmentMetricVolume, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	alignmentMetricVolume = flip(alignmentMetricVolume,1);
end
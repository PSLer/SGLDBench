function alignmentMetricVolume = Common_ComputeEdgeAlignmentDeviation(dominantDirDesign)
	global meshHierarchy_;
	global densityLayout_;
	global voxelsOnBoundary_;
	global frameStruct4Voxelization_; 
	
	numElements = meshHierarchy_(1).numElements;
	alignmentMetric = zeros(size(dominantDirDesign,1),1);
	
	solidElementsInDesign = densityLayout_>=0.1;
	assoVoxels = frameStruct4Voxelization_.AssociatedVoxels;
	graphNodeCoords = frameStruct4Voxelization_.nodeCoords;
	graphEdges = frameStruct4Voxelization_.eNodMat;
	
	for ii=1:numel(assoVoxels)
		iEdge = assoVoxels(ii).iEdge;
		iEdgeEndPots = graphNodeCoords(graphEdges(iEdge,:),:);
		iEdgeDir = iEdgeEndPots(1,:) - iEdgeEndPots(2,:);
		iEdgeDir = iEdgeDir / norm(iEdgeDir);
		iStressDomiDirs = dominantDirDesign(assoVoxels(ii).arr,:);
		iEdgeACOS = iStressDomiDirs * iEdgeDir(:) ./ vecnorm(iStressDomiDirs,2,2);
		iEdgeAngSin = sqrt(1 - iEdgeACOS.^2);
		alignmentMetric(assoVoxels(ii).arr,1) = iEdgeAngSin;
alignmentMetric(assoVoxels(ii).arr,1) = 1 - iEdgeAngSin;		
	end
	alignmentMetric(voxelsOnBoundary_) = 0;
	alignmentMetricVolume = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	alignmentMetricVolume(meshHierarchy_(1).eleMapBack,1) = alignmentMetric;
	alignmentMetricVolume = reshape(alignmentMetricVolume, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	alignmentMetricVolume = flip(alignmentMetricVolume,1);
end
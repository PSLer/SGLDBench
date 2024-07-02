function volumeFractionOfVoxelizedMeshEdges = MGD_CoatMeshEdgesWithVoxels(edgeWidth, numLayerboundary, numLayerLoads, numLayerFixation)
	global meshHierarchy_;
	global voxelizedMeshEdges_;
	global voxelizedMeshEdgesWithGivenWidth_;
	global densityLayoutWithoutBoundary_;
	global densityLayout_;
	
	%%Mesh Edges
	voxelizedMeshEdgesWithGivenWidth_ = voxelizedMeshEdges_;
	densityLayoutWithoutBoundary_ = zeros(meshHierarchy_(1).numElements,1);
	for ii=1:edgeWidth
		blockIndex = MissionPartition(numel(voxelizedMeshEdgesWithGivenWidth_), 1.0e7);
		numBlocks = size(blockIndex,1);
		iVoxels = struct('arr', []); iVoxels = repmat(iVoxels, numBlocks);
		for jj=1:numBlocks
			iVoxels(jj).arr = voxelizedMeshEdgesWithGivenWidth_(blockIndex(jj,1):blockIndex(jj,2),:)';
			iVoxels(jj).arr = Common_IncludeAdjacentElements(iVoxels(jj).arr);
		end
		voxelizedMeshEdgesWithGivenWidth_ = [iVoxels.arr]; 
        voxelizedMeshEdgesWithGivenWidth_ = voxelizedMeshEdgesWithGivenWidth_(:);
	end	
	[boundaryElementsWithGivenThickness, passiveElementsNearLoads, passiveElementsNearFixation] = ...
		TopOpti_SetPassiveElements(numLayerboundary, numLayerLoads, numLayerFixation);
	voxelizedMeshEdgesWithGivenWidth_ = unique([voxelizedMeshEdgesWithGivenWidth_; passiveElementsNearLoads; passiveElementsNearFixation]);
	densityLayoutWithoutBoundary_(voxelizedMeshEdgesWithGivenWidth_,1) = 1;
	voxelizedMeshEdgesWithGivenWidth_ = unique([voxelizedMeshEdgesWithGivenWidth_; boundaryElementsWithGivenThickness]);
	densityLayout_(voxelizedMeshEdgesWithGivenWidth_,1) = 1;

	volumeFractionOfVoxelizedMeshEdges = numel(voxelizedMeshEdgesWithGivenWidth_) / meshHierarchy_(1).numElements;
	% disp(['Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionOfVoxelizedMeshEdges)]);
	% totVolume = meshHierarchy_(1).numElements * (meshHierarchy_(1).eleSize(1)*meshHierarchy_(1).eleSize(2)*meshHierarchy_(1).eleSize(3));
    % frameVolume = numel(voxelizedMeshEdgesWithGivenWidth_) * (meshHierarchy_(1).eleSize(1)*meshHierarchy_(1).eleSize(2)*meshHierarchy_(1).eleSize(3));
    % disp(['Frame Volume of Mesh Edges: ', sprintf('%.6f', frameVolume), ' of ', sprintf('%.6f', totVolume)]);	
end

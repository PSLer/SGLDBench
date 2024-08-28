function voxelsAlongLatticeEdges = MGD_CoatMeshEdgesWithVoxels_B(edgeWidth, passiveElements)
	global meshHierarchy_;
	global voxelizedMeshEdges_;
	global densityLayout_;
	
	%%Mesh Edges
	voxelsAlongLatticeEdges = voxelizedMeshEdges_;
	for ii=1:edgeWidth-1
		blockIndex = Solving_MissionPartition(numel(voxelsAlongLatticeEdges), 1.0e7);
		numBlocks = size(blockIndex,1);
		iVoxels = struct('arr', []); iVoxels = repmat(iVoxels, numBlocks);
		for jj=1:numBlocks
			iVoxels(jj).arr = voxelsAlongLatticeEdges(blockIndex(jj,1):blockIndex(jj,2),:)';
			iVoxels(jj).arr = Common_IncludeAdjacentElements(iVoxels(jj).arr);
            iVoxels(jj).arr = iVoxels(jj).arr(:)';
		end
		voxelsAlongLatticeEdges = [iVoxels.arr]; 
        voxelsAlongLatticeEdges = voxelsAlongLatticeEdges(:);
	end	

	voxelsAlongLatticeEdges = unique([voxelsAlongLatticeEdges(:); passiveElements(:)]);
	% densityLayout_(voxelsAlongLatticeEdges,1) = 1;

	% volumeFractionOfVoxelizedMeshEdges = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
	% disp(['Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionOfVoxelizedMeshEdges)]);
	% totVolume = meshHierarchy_(1).numElements * (meshHierarchy_(1).eleSize(1)*meshHierarchy_(1).eleSize(2)*meshHierarchy_(1).eleSize(3));
    % frameVolume = numel(voxelsAlongLatticeEdges) * (meshHierarchy_(1).eleSize(1)*meshHierarchy_(1).eleSize(2)*meshHierarchy_(1).eleSize(3));
    % disp(['Frame Volume of Mesh Edges: ', sprintf('%.6f', frameVolume), ' of ', sprintf('%.6f', totVolume)]);	
end

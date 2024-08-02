function volumeFractionOfVoxelizedMeshEdges = PSLs_VoxelizePSLs(numLayerPSLs, numLayerboundary, numLayerLoads, numLayerFixation)
	global meshHierarchy_;
	global majorPSLpool_;
	global mediumPSLpool_;
	global minorPSLpool_;
	global loadingCond_;
	global fixingCond_;
	
	global densityLayout_;
	global densityLayoutWithoutBoundary_;
	global PSLs2Bvoxelized_;
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;
	
	densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	miniPSLength = 20;
	tarMajorPSLs = majorPSLpool_;
	tarIndice = [];
	for ii=1:length(majorPSLpool_)
		if tarMajorPSLs(ii).length > miniPSLength
			tarIndice(end+1,1) = ii;
			tarMajorPSLs(ii).eleIndexList = tarMajorPSLs(ii).eleIndexList(:)';
		end
	end
	tarMajorPSLs = tarMajorPSLs(tarIndice);
	
	tarMediumPSLs = mediumPSLpool_;
	tarIndice = [];
	for ii=1:length(tarMediumPSLs)
		if tarMediumPSLs(ii).length > miniPSLength
			tarIndice(end+1,1) = ii;
			tarMediumPSLs(ii).eleIndexList = tarMediumPSLs(ii).eleIndexList(:)';
		end
	end
	tarMediumPSLs = tarMediumPSLs(tarIndice);	

	tarMinorPSLs = minorPSLpool_;
	tarIndice = [];
	for ii=1:length(tarMinorPSLs)
		if tarMinorPSLs(ii).length > miniPSLength
			tarIndice(end+1,1) = ii;
			tarMinorPSLs(ii).eleIndexList = tarMinorPSLs(ii).eleIndexList(:)';
		end
	end
	tarMinorPSLs = tarMinorPSLs(tarIndice);

	PSLs2Bvoxelized_ = [tarMajorPSLs; tarMediumPSLs; tarMinorPSLs];
	
	%%On PSLs
	voxelsAlongPSLs = [PSLs2Bvoxelized_.eleIndexList]; voxelsAlongPSLs = voxelsAlongPSLs(:);
	for ii=1:numLayerPSLs-1
		blockIndex = Solving_MissionPartition(numel(voxelsAlongPSLs), 1.0e7);
		numBlocks = size(blockIndex,1);
		iVoxels = struct('arr', []); iVoxels = repmat(iVoxels, numBlocks);
		for jj=1:numBlocks
			iVoxels(jj).arr = voxelsAlongPSLs(blockIndex(jj,1):blockIndex(jj,2),:)';
			iVoxels(jj).arr = Common_IncludeAdjacentElements(iVoxels(jj).arr);
		end
		voxelsAlongPSLs = [iVoxels.arr]'; voxelsAlongPSLs = voxelsAlongPSLs(:);
	end
	
	[voxelsOnBoundary_, voxelsInLoadingArea_, voxelsInFixingArea_] = TopOpti_SetPassiveElements(numLayerboundary, numLayerLoads, numLayerFixation);
	densityLayout_(voxelsAlongPSLs,1) = 1;
	densityLayout_(voxelsInLoadingArea_,1) = 1;
	densityLayout_(voxelsInFixingArea_,1) = 1;
	densityLayoutWithoutBoundary_ = densityLayout_;
	densityLayout_(voxelsOnBoundary_,1) = 1;
	
	volumeFractionOfVoxelizedMeshEdges = sum(densityLayout_) / meshHierarchy_(1).numElements;
	disp(['Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionOfVoxelizedMeshEdges)]);	
end
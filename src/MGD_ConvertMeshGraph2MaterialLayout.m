function MGD_ConvertMeshGraph2MaterialLayout(edgeWidth, numLayerboundary, numLayerLoads, numLayerFixation)
	global meshHierarchy_;
	global volumeFractionDesign_; 
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;	
	global densityLayout_;
	global densityLayout4Vis_;
	

	densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	densityLayout4Vis_ = zeros(size(meshHierarchy_(1).eleMapForward));
	
	tStart = tic;
	[voxelsOnBoundary_, voxelsInLoadingArea_, voxelsInFixingArea_] = TopOpti_SetPassiveElements(numLayerboundary, numLayerLoads, numLayerFixation);
	passiveElements = unique([voxelsOnBoundary_(:); voxelsInLoadingArea_(:); voxelsInFixingArea_(:)]);	
	
	[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);
	volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
	disp(sprintf('Design Volume Fraction %.6f', volumeFractionDesign_));			
	densityLayout_(voxelsAlongLatticeEdges) = 1; 
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
	densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;		
	
	tEnd = toc(tStart);
	disp(['............Convert Mesh/Graph to Material Layout Costs: ', sprintf('%.1f', tEnd), 's']);	
end

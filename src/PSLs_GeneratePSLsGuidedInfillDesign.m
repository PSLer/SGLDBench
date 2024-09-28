function PSLs_GeneratePSLsGuidedInfillDesign(psDirIndicator, edgeWidth, targetDepositionRatio, numLayerboundary, numLayerLoads, numLayerFixation)
	global meshHierarchy_;
	global volumeFractionDesign_;
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;
	global densityLayout_;
	global densityLayout4Vis_;
	
	upperLineDensCtrl = 30;
	lowerLineDensCtrl = 10;
	permittedVolumeDeviation = 0.03;
	opt_DetermingUpperBound = 1;
	densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	densityLayout4Vis_ = zeros(size(meshHierarchy_(1).eleMapForward));
	if targetDepositionRatio>0.9
		warning('Close to a solid domain, no need for design!');
		densityLayout_ = ones(size(densityLayout_));
		densityLayout4Vis_(meshHierarchy_(1).eleMapBack) = 1;
		volumeFractionDesign_ = 1;
		tEnd = toc(tStart);
		disp(['Conduct PSLs-guided Structural Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);				
		return;
	end
	
	tStart = tic;
	[voxelsOnBoundary_, voxelsInLoadingArea_, voxelsInFixingArea_] = TopOpti_SetPassiveElements(numLayerboundary, numLayerLoads, numLayerFixation);
	passiveElements = unique([voxelsOnBoundary_(:); voxelsInLoadingArea_(:); voxelsInFixingArea_(:)]);
	
	%% Check Design Space
	volumeFractionDesign_ = numel(passiveElements) / meshHierarchy_(1).numElements;
	if volumeFractionDesign_ > targetDepositionRatio
		disp(['Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionDesign_)]);	
		warning('Too many passive elements, there is no design space!');
		densityLayout_(passiveElements,1) = 1;
		densityLayout4Vis_(meshHierarchy_(1).eleMapBack(passiveElements),1) = 1;
		tEnd = toc(tStart);
		disp(['Conduct PSLs-guided Structural Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);				
		return;
	end
	
	%% Determine the lower bound for PSL density control
	volumeFractionDesign_ = 1;
	while volumeFractionDesign_ > targetDepositionRatio
		lineDensCtrl = lowerLineDensCtrl;
		PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator);
		PSLs_ConvertPSLs2PiecewiseGraphs();
		[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);						
		voxelsOutOfVolume = setdiff(voxelsAlongLatticeEdgesWithoutPassiveElesMapback, meshHierarchy_(1).eleMapBack);
		voxelsAlongLatticeEdgesWithoutPassiveElesMapback = setdiff(voxelsAlongLatticeEdgesWithoutPassiveElesMapback, voxelsOutOfVolume);
		
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['Determining Lower Bound of PSL Density Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Line Density Para %.1f', lineDensCtrl)]);
		if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
			densityLayout_(voxelsAlongLatticeEdges) = 1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
			densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;			
			tEnd = toc(tStart);
			disp(['Conduct PSLs-guided Structural Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);					
			return;
		end			
		if volumeFractionDesign_ > targetDepositionRatio
			upperLineDensCtrl = lowerLineDensCtrl; opt_DetermingUpperBound = 0;
			lowerLineDensCtrl = lowerLineDensCtrl / 1.25;
		else
			break;
		end
		if lowerLineDensCtrl < 2
			warning('Inappropriate settings for the material budget!');
			densityLayout_(voxelsAlongLatticeEdges) = 1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
			densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;
			tEnd = toc(tStart);
			disp(['Conduct PSLs-guided Structural Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);					
			return;
		end
	end
	
	%% Determine the upper bound for PSL density control
	if opt_DetermingUpperBound
		volumeFractionDesign_ = 0;
		while volumeFractionDesign_ < targetDepositionRatio
			lineDensCtrl = upperLineDensCtrl;
			PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator);
			PSLs_ConvertPSLs2PiecewiseGraphs();
			[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);						
			voxelsOutOfVolume = setdiff(voxelsAlongLatticeEdgesWithoutPassiveElesMapback, meshHierarchy_(1).eleMapBack);
			voxelsAlongLatticeEdgesWithoutPassiveElesMapback = setdiff(voxelsAlongLatticeEdgesWithoutPassiveElesMapback, voxelsOutOfVolume);
			
			volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
			disp(['Determining Upper Bound of PSL Density Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
				sprintf(' with Line Density Para %.1f', lineDensCtrl)]);
			if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
				densityLayout_(voxelsAlongLatticeEdges) = 1; 
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
				densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;				
				tEnd = toc(tStart);
				disp(['Conduct PSLs-guided Structural Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);				
				return;
			end			
			if volumeFractionDesign_ < targetDepositionRatio
				lowerLineDensCtrl = upperLineDensCtrl;
				upperLineDensCtrl = upperLineDensCtrl * 1.25;
			else
				break;
			end		
		end
	end	

	
	%%Determine the target PSL density control
	idx = 1;
	while abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio > permittedVolumeDeviation			
		lineDensCtrl = (lowerLineDensCtrl + upperLineDensCtrl) / 2;
		PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator);
		PSLs_ConvertPSLs2PiecewiseGraphs();
		[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);						
		voxelsOutOfVolume = setdiff(voxelsAlongLatticeEdgesWithoutPassiveElesMapback, meshHierarchy_(1).eleMapBack);
		voxelsAlongLatticeEdgesWithoutPassiveElesMapback = setdiff(voxelsAlongLatticeEdgesWithoutPassiveElesMapback, voxelsOutOfVolume);
			
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['Design Iteration ', sprintf('%d', idx), sprintf('. Design Volume Fraction: %.6f', volumeFractionDesign_), ...
			sprintf(' with Line Density Para %.1f', lineDensCtrl)]);
		if volumeFractionDesign_>targetDepositionRatio
			upperLineDensCtrl = lineDensCtrl;
		else
			lowerLineDensCtrl = lineDensCtrl;
		end
		idx = idx + 1;
		if idx > 10
			warning('PSLs-guided Infill failed to converge to the prescribed design'); break;
		end
	end	
	densityLayout_(voxelsAlongLatticeEdges) = 1;
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
	densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;
	tEnd = toc(tStart);
	disp(['Conduct PSLs-guided Structural Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);
end

function voxelsAlongLatticeEdges = PSLs_GetVoxelsPassedByPSLs(edgeWidth, passiveElements)
	global meshHierarchy_;
	global majorPSLpool_;
	global mediumPSLpool_;
	global minorPSLpool_;
	global loadingCond_;
	global fixingCond_;
	
	global densityLayout_;
	global PSLs2Bvoxelized_;
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;
	
	% densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
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
	voxelsAlongLatticeEdges = [PSLs2Bvoxelized_.eleIndexList]; voxelsAlongLatticeEdges = unique(voxelsAlongLatticeEdges(:));
	for ii=1:edgeWidth-1
		blockIndex = Solving_MissionPartition(numel(voxelsAlongLatticeEdges), 1.0e7);
		numBlocks = size(blockIndex,1);
		iVoxels = struct('arr', []); iVoxels = repmat(iVoxels, numBlocks);
		for jj=1:numBlocks
			iVoxels(jj).arr = voxelsAlongLatticeEdges(blockIndex(jj,1):blockIndex(jj,2),:)';
			iVoxels(jj).arr = Common_IncludeAdjacentElements(iVoxels(jj).arr);
            iVoxels(jj).arr = iVoxels(jj).arr(:)';
		end
		voxelsAlongLatticeEdges = [iVoxels.arr]'; 
        voxelsAlongLatticeEdges = voxelsAlongLatticeEdges(:);
	end
	voxelsAlongLatticeEdges = unique([voxelsAlongLatticeEdges(:); passiveElements(:)]);
	% densityLayout_(voxelsAlongLatticeEdges,1) = 1;
	% densityLayout_(passiveElements,1) = 1
	
	% volumeFractionOfVoxelizedMeshEdges = sum(densityLayout_) / meshHierarchy_(1).numElements;
	% disp(['Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionOfVoxelizedMeshEdges)]);	
end

function PSLs_ConvertPSLs2PiecewiseGraphs(piecewiseSpan)
	global majorPSLpool_;
	global mediumPSLpool_;
	global minorPSLpool_;
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	frameStruct4Voxelization_ = Data_VertexEdgeGraphStruct();

	% densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	miniPSLength = 20;
	graphNodeCoords = [];
	graphEdges = [];
	for ii=1:numel(majorPSLpool_)
		if majorPSLpool_(ii).length > miniPSLength
			numExisingNode = size(graphNodeCoords,1);
			[iGraphNodes, iGraphEdges] = PSLs_PartitionPSL2Graph(majorPSLpool_(ii).phyCoordList, miniPSLength, numExisingNode);
			graphNodeCoords(end+1:end+size(iGraphNodes,1),:) = iGraphNodes;
			graphEdges(end+1:end+size(iGraphEdges,1),:) = iGraphEdges;
		end	
	end
	for ii=1:numel(mediumPSLpool_)
		if mediumPSLpool_(ii).length > miniPSLength
			numExisingNode = size(graphNodeCoords,1);
			[iGraphNodes, iGraphEdges] = PSLs_PartitionPSL2Graph(mediumPSLpool_(ii).phyCoordList, miniPSLength, numExisingNode);
			graphNodeCoords(end+1:end+size(iGraphNodes,1),:) = iGraphNodes;
			graphEdges(end+1:end+size(iGraphEdges,1),:) = iGraphEdges;
		end	
	end		
	for ii=1:numel(minorPSLpool_)
		if minorPSLpool_(ii).length > miniPSLength
			numExisingNode = size(graphNodeCoords,1);
			[iGraphNodes, iGraphEdges] = PSLs_PartitionPSL2Graph(minorPSLpool_(ii).phyCoordList, miniPSLength, numExisingNode);
			graphNodeCoords(end+1:end+size(iGraphNodes,1),:) = iGraphNodes;
			graphEdges(end+1:end+size(iGraphEdges,1),:) = iGraphEdges;
		end	
	end
	
	frameStruct4Voxelization_.numNodes = size(graphNodeCoords,1);
	frameStruct4Voxelization_.nodeCoords = graphNodeCoords;
	frameStruct4Voxelization_.numEdges = size(graphEdges,1);
	frameStruct4Voxelization_.eNodMat = graphEdges;
	frameStruct4Voxelization_.edgeLengths = vecnorm(frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,1),:) ...
		- frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,2),:),2,2);	
	vertexEdgeGraph_ = frameStruct4Voxelization_;
end


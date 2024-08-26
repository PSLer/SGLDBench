function PSLs_GeneratePSLsGuidedInfillDesign(psDirIndicator, numLayerPSLs, targetDepositionRatio, numLayerboundary, numLayerLoads, numLayerFixation)
	global boundingBox_;
	global meshHierarchy_;
	global volumeFractionDesign_;
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;
	global densityLayout_;

	upperLineDensCtrl = 20;
	lowerLineDensCtrl = 5;
	permittedVolumeDeviation = 0.05;
	
	densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	if targetDepositionRatio>0.9
		warning('Close to a solid domain, no need for design!');
		densityLayout_ = ones(size(densityLayout_));
		volumeFractionDesign_ = 1;
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
		return;
	end
	
	%% Determine the lower bound for PSL density control
	volumeFractionDesign_ = 1;
	while volumeFractionDesign_ > targetDepositionRatio
		lineDensCtrl = lowerLineDensCtrl;
		PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator);
		voxelsAlongPSLs = PSLs_GetVoxelsPassedByPSLs(numLayerPSLs, passiveElements);
		volumeFractionDesign_ = numel(voxelsAlongPSLs) / meshHierarchy_(1).numElements;
		disp(['Determining Lower Bound of PSL Density Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Line Density Para %.1f', lineDensCtrl)]);
		if volumeFractionDesign_ > targetDepositionRatio
			lowerLineDensCtrl = lowerLineDensCtrl / 1.25;
		else
			break;
		end
		if lowerLineDensCtrl < 2
			warning('Inappropriate settings for the material budget!');
			densityLayout_(voxelsAlongPSLs) = 1;
			return;
		end
	end
	if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
		densityLayout_(voxelsAlongPSLs) = 1; return;
	end
	
	%% Determine the upper bound for PSL density control
	volumeFractionDesign_ = 0;
	while volumeFractionDesign_ < targetDepositionRatio
		lineDensCtrl = upperLineDensCtrl;
		PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator);
		voxelsAlongPSLs = PSLs_GetVoxelsPassedByPSLs(numLayerPSLs, passiveElements);
		volumeFractionDesign_ = numel(voxelsAlongPSLs) / meshHierarchy_(1).numElements;
		disp(['Determining Upper Bound of PSL Density Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Line Density Para %.1f', lineDensCtrl)]);
		if volumeFractionDesign_ < targetDepositionRatio
			upperLineDensCtrl = upperLineDensCtrl * 1.25;
		else
			break;
		end		
	end
	if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
		densityLayout_(voxelsAlongPSLs) = 1; return;
	end
	
	%%Determine the target PSL density control
	idx = 1;
	while abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio > permittedVolumeDeviation			
		lineDensCtrl = (lowerLineDensCtrl + upperLineDensCtrl) / 2;
		PSLs_GeneratePSLsBy3DTSV(lineDensCtrl, psDirIndicator);
		voxelsAlongPSLs = PSLs_GetVoxelsPassedByPSLs(numLayerPSLs, passiveElements);
		volumeFractionDesign_ = numel(voxelsAlongPSLs) / meshHierarchy_(1).numElements;
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
	densityLayout_(voxelsAlongPSLs) = 1;
	tEnd = toc(tStart);
	disp(['Conduct PSLs-guided Structural Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);
end

function voxelsAlongPSLs = PSLs_GetVoxelsPassedByPSLs(numLayerPSLs, passiveElements)
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
	voxelsAlongPSLs = [PSLs2Bvoxelized_.eleIndexList]; voxelsAlongPSLs = voxelsAlongPSLs(:);
	for ii=1:numLayerPSLs-1
		blockIndex = Solving_MissionPartition(numel(voxelsAlongPSLs), 1.0e7);
		numBlocks = size(blockIndex,1);
		iVoxels = struct('arr', []); iVoxels = repmat(iVoxels, numBlocks);
		for jj=1:numBlocks
			iVoxels(jj).arr = voxelsAlongPSLs(blockIndex(jj,1):blockIndex(jj,2),:)';
			iVoxels(jj).arr = Common_IncludeAdjacentElements(iVoxels(jj).arr);
            iVoxels(jj).arr = iVoxels(jj).arr(:)';
		end
		voxelsAlongPSLs = [iVoxels.arr]'; 
        voxelsAlongPSLs = voxelsAlongPSLs(:);
	end
	voxelsAlongPSLs = unique([voxelsAlongPSLs(:); passiveElements(:)]);
	% densityLayout_(voxelsAlongPSLs,1) = 1;
	% densityLayout_(passiveElements,1) = 1
	
	% volumeFractionOfVoxelizedMeshEdges = sum(densityLayout_) / meshHierarchy_(1).numElements;
	% disp(['Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionOfVoxelizedMeshEdges)]);	
end
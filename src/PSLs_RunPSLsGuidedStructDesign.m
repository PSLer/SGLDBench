function volumeFractionOfVoxelizedMeshEdges = PSLs_RunPSLsGuidedStructDesign(targetDepositionRatio, PSLthickness)
	global meshHierarchy_;
	global strainEnergyPerElement_;
	
	global PSLs2Bvoxelized_;
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;
	global densityLayout_;
	global densityLayoutWithoutBoundary_;
	
	densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	%%1. setup
	%%1.1 elements along PSLs
	numTrajectories = numel(PSLs2Bvoxelized_);
	for ii=1:numTrajectories
		iPSL = PSLs2Bvoxelized_(ii);		
		elesAroundPSL = iPSL.eleIndexList;
		index = 2;
		while index<=PSLthickness
			elesAroundPSL = Common_IncludeAdjacentElements(elesAroundPSL);
			index = index + 1;
		end
		PSLs2Bvoxelized_(ii).adjacentVoxels = elesAroundPSL(:)';
	end
	opt = 1;
	voxelsAlongPSLs = [PSLs2Bvoxelized_.adjacentVoxels];
	
	numTargetElements = round(meshHierarchy_(1).numElements*targetDepositionRatio);
	solidVoxels = [voxelsAlongPSLs(:); voxelsInLoadingArea_(:); voxelsInFixingArea_(:); voxelsOnBoundary_(:)];
	solidVoxels = unique(solidVoxels);
	numRealElements = numel(solidVoxels);	
	
	if solidVoxels >= numTargetElements
		opt = 0;
		densityLayout_(voxelsAlongPSLs,1) = 1;	
		densityLayout_(voxelsInLoadingArea_,1) = 1;
		densityLayout_(voxelsInFixingArea_,1) = 1;
		densityLayoutWithoutBoundary_ = densityLayout_;
		densityLayout_(voxelsOnBoundary_,1) = 1;		
		warning('Insufficient Space for Postprocess!');
		return;
	end
	
	%%1.3 Approximate metric for compliance
	for ii=1:numTrajectories
		PSLs2Bvoxelized_(ii).strainEnergyApprox = sum(strainEnergyPerElement_(PSLs2Bvoxelized_(ii).adjacentVoxels));
	end
	complianceBasedImportanceMetric = [PSLs2Bvoxelized_.strainEnergyApprox];
	complianceBasedImportanceMetric = complianceBasedImportanceMetric(:)';
	
	%%2. find skeleton PSL bands with proper thickness value 
	basisElementBands = struct('coreElements', [], 'coatingElements', []);
	basisElementBands = repmat(basisElementBands, numTrajectories, 1);
	for ii=1:numTrajectories
		basisElementBands(ii).coreElements = PSLs2Bvoxelized_(ii).adjacentVoxels;
	end
	while numRealElements<numTargetElements
		for ii=1:numTrajectories
			iNewCoreElements = [basisElementBands(ii).coreElements basisElementBands(ii).coatingElements];
			iNewBand = Common_IncludeAdjacentElements(iNewCoreElements);
			iNewBand = iNewBand(:)';
			% if strcmp(paraType, 'PSL Lattice Structure')
				% iNewBand = intersect(iNewBand, PSLs2Bvoxelized_(ii).accompEles);
			% end
			iNewCoatingElements = setdiff(iNewBand, iNewCoreElements);
			basisElementBands(ii).coreElements = iNewCoreElements;
			basisElementBands(ii).coatingElements = iNewCoatingElements;
		end
		solidVoxels = [[basisElementBands.coreElements], [basisElementBands.coatingElements], ...
			voxelsOnBoundary_(:)', voxelsInLoadingArea_(:)', voxelsInFixingArea_(:)'];
		solidVoxels = unique(solidVoxels);
		numRealElements = numel(solidVoxels);
	end
	thicknessUniformityOpt = 0;
	if ~thicknessUniformityOpt
		for ii=1:numTrajectories
			iAllElements = unique([basisElementBands(ii).coreElements basisElementBands(ii).coatingElements]);
			basisElementBands(ii).coreElements = PSLs2Bvoxelized_(ii).adjacentVoxels;
			iAllElementsAdded = Common_IncludeAdjacentElements(iAllElements); 
			iAllElements = iAllElementsAdded(:)';
			% if strcmp(paraType, 'PSL Lattice Structure')
				% iAllElements = intersect(iAllElements, PSLs2Bvoxelized_(ii).accompEles);
			% end			
			basisElementBands(ii).coatingElements = setdiff(iAllElements,basisElementBands(ii).coreElements);
		end
	end
	
	%%3. Adjust thickness per trajectory
	[~, reOrderMap] = sort(complianceBasedImportanceMetric, 'descend');
	basisElementBands = basisElementBands(reOrderMap);
	PSLs2Bvoxelized_ = PSLs2Bvoxelized_(reOrderMap);
	for ii=1:numTrajectories
		PSLs2Bvoxelized_(ii).adjacentVoxels = basisElementBands(ii).coreElements;
	end
	solidVoxels = unique([[basisElementBands.coreElements], voxelsOnBoundary_(:)', voxelsInLoadingArea_(:)', voxelsInFixingArea_(:)']);
	for ii=1:numTrajectories	
		tmpSolidVoxels = unique([solidVoxels basisElementBands(ii).coatingElements]);
		if numel(tmpSolidVoxels)>=numTargetElements, break; end
		solidVoxels = tmpSolidVoxels;
		PSLs2Bvoxelized_(ii).adjacentVoxels = unique([basisElementBands(ii).coreElements basisElementBands(ii).coatingElements]);
	end
	
	densityLayout_(solidVoxels,1) = 1;
	densityLayoutWithoutBoundary_ = densityLayout_;
	densityLayoutWithoutBoundary_(voxelsOnBoundary_,1) = 0;
	densityLayoutWithoutBoundary_(voxelsInLoadingArea_,1) = 1;
	densityLayoutWithoutBoundary_(voxelsInFixingArea_,1) = 1;
	
	volumeFractionOfVoxelizedMeshEdges = sum(densityLayout_) / meshHierarchy_(1).numElements;
	disp(['Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionOfVoxelizedMeshEdges)]);		
end


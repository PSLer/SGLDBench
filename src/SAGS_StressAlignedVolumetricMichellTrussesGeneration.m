function SAGS_StressAlignedVolumetricMichellTrussesGeneration(edgeWidth, targetDepositionRatio, numLayerboundary, ...
                numLayerLoads, numLayerFixation)
	global meshHierarchy_;
	global volumeFractionDesign_; 
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;	
	global densityLayout_;
	global optEdgeAlignmentComparison_; optEdgeAlignmentComparison_ = 1;
	
	upperLatticeSizeCtrl = 32;
	lowerLatticeSizeCtrl = 12;
	
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
		disp(['............Volume Fraction of Mesh Edges: ' sprintf('%16.6g',volumeFractionDesign_)]);	
		warning('Too many passive elements, there is no design space!');
		densityLayout_(passiveElements,1) = 1;
		return;
	end	
	
	%% Determine the upper bound for lattice size control
	volumeFractionDesign_ = 0;
	while volumeFractionDesign_ < targetDepositionRatio
		latticeSizeCtrl = round(upperLatticeSizeCtrl);

		SAGS_CallArora2019MatlabSuite_ExtractingGraph(latticeSizeCtrl);

		%Voxelization%
		if optEdgeAlignmentComparison_
			voxelsAlongLatticeEdges = MGD_VoxelizeMeshEdges_PerEdge(edgeWidth, passiveElements);		
		else
			MGD_VoxelizeMeshEdges();
			voxelsAlongLatticeEdges = MGD_CoatMeshEdgesWithVoxels_B(edgeWidth, passiveElements);				
		end
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Determining Upper Bound of Lattice Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Size Ctrl Para %.1f', latticeSizeCtrl)]);
		if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
			densityLayout_(voxelsAlongLatticeEdges) = 1; return;
		end		
		if volumeFractionDesign_ < targetDepositionRatio
			upperLatticeSizeCtrl = upperLatticeSizeCtrl * 1.25;
		else
			break;
		end
		if upperLatticeSizeCtrl > 64
			warning('Inappropriate settings for the material budget!');
			densityLayout_(voxelsAlongLatticeEdges) = 1;
			return;
		end		
	end	
	
	%% Determine the lower bound for lattice size control
	volumeFractionDesign_ = 1;
	while volumeFractionDesign_ > targetDepositionRatio
		latticeSizeCtrl = round(lowerLatticeSizeCtrl);
		
		SAGS_CallArora2019MatlabSuite_ExtractingGraph(latticeSizeCtrl);

		%Voxelization%
		if optEdgeAlignmentComparison_
			voxelsAlongLatticeEdges = MGD_VoxelizeMeshEdges_PerEdge(edgeWidth, passiveElements);		
		else
			MGD_VoxelizeMeshEdges();
			voxelsAlongLatticeEdges = MGD_CoatMeshEdgesWithVoxels_B(edgeWidth, passiveElements);				
		end
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Determining Lower Bound of Lattice Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Size Ctrl Para %.1f', latticeSizeCtrl)]);		
		if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
			densityLayout_(voxelsAlongLatticeEdges) = 1; return;
		end			
		if volumeFractionDesign_ > targetDepositionRatio
			lowerLatticeSizeCtrl = lowerLatticeSizeCtrl / 1.25;
		else
			break;
		end		
		if lowerLatticeSizeCtrl < 6
			warning('Inappropriate settings for the material budget!');
			densityLayout_(voxelsAlongLatticeEdges) = 1;
			return;
		end			
	end	

	%%Determine the target Lattice Size density control
	idx = 1;
	while abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio > permittedVolumeDeviation			
		latticeSizeCtrl = round((lowerLatticeSizeCtrl + upperLatticeSizeCtrl) / 2);
		
		SAGS_CallArora2019MatlabSuite_ExtractingGraph(latticeSizeCtrl);

		%Voxelization%
		if optEdgeAlignmentComparison_
			voxelsAlongLatticeEdges = MGD_VoxelizeMeshEdges_PerEdge(edgeWidth, passiveElements);		
		else
			MGD_VoxelizeMeshEdges();
			voxelsAlongLatticeEdges = MGD_CoatMeshEdgesWithVoxels_B(edgeWidth, passiveElements);				
		end
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Design Iteration ', sprintf('%d', idx), sprintf('. Design Volume Fraction: %.6f', volumeFractionDesign_), ...
			sprintf(' with Line Density Para %.1f', latticeSizeCtrl)]);
		
		if volumeFractionDesign_>targetDepositionRatio					
            upperLatticeSizeCtrl = latticeSizeCtrl;
		else
			lowerLatticeSizeCtrl = latticeSizeCtrl;
		end
		idx = idx + 1;
		if idx > 10
			warning('Stress-aligned Volumetric Michell Truss Infill failed to converge to the prescribed design'); break;
		end
	end	
	densityLayout_(voxelsAlongLatticeEdges) = 1;
	tEnd = toc(tStart);
	disp(['............Conduct Stress-aligned Volumetric Michell Trusses Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);	
end

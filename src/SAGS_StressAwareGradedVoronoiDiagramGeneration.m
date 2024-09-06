function SAGS_StressAwareGradedVoronoiDiagramGeneration(edgeWidth, targetDepositionRatio, numLayerboundary, ...
                numLayerLoads, numLayerFixation, sizeAspectRatio)
	global meshHierarchy_;
	global volumeFractionDesign_; 
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;	
	global densityLayout_;
	global optEdgeAlignmentComparison_; optEdgeAlignmentComparison_ = 1;
	
	upperLatticeSizeCtrl = 0.1;
	lowerLatticeSizeCtrl = 0.04;
	
	permittedVolumeDeviation = 0.05;
	opt_DetermingLowerBound = 1;
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
	volumeFractionDesign_ = 1;
	while volumeFractionDesign_ > targetDepositionRatio
		latticeSizeCtrl = upperLatticeSizeCtrl;
	
		%%Create Graded Voronoi Diagram
		callGradedVoronoiGenerater_python = ['./externalModules/GradedVoronoiDiagram/AdaptiveGraphGenerator.py ./out/StressField_Tet_v2.stress ', ...
			sprintf('%g ',latticeSizeCtrl), sprintf('%g',sizeAspectRatio)];
		pyrunfile(callGradedVoronoiGenerater_python);
		LoadGeneratedGraphFromFileObj_B();
			
		%Voxelization%
		if optEdgeAlignmentComparison_
			[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveEles] = MGD_VoxelizeMeshEdges_PerEdge(edgeWidth, passiveElements);		
		else
			MGD_VoxelizeMeshEdges();
			voxelsAlongLatticeEdges = MGD_CoatMeshEdgesWithVoxels_B(edgeWidth, passiveElements);				
		end
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Determining Upper Bound of Lattice Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Size Ctrl Para %g', latticeSizeCtrl)]);
		if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
			densityLayout_(voxelsAlongLatticeEdges) = 1; return;
		end		
		if volumeFractionDesign_ > targetDepositionRatio
			lowerLatticeSizeCtrl = upperLatticeSizeCtrl;
			opt_DetermingLowerBound = 0;
			upperLatticeSizeCtrl = upperLatticeSizeCtrl * 1.25;
		else
			break;
		end
		if upperLatticeSizeCtrl > 0.5
			warning('Inappropriate settings for the material budget!');
			densityLayout_(voxelsAlongLatticeEdges) = 1;
			return;
		end		
	end

	if opt_DetermingLowerBound
		%% Determine the lower bound for lattice size control
		volumeFractionDesign_ = 0;
		while volumeFractionDesign_ < targetDepositionRatio
			latticeSizeCtrl = lowerLatticeSizeCtrl;
	
			%%Create Graded Voronoi Diagram
			callGradedVoronoiGenerater_python = ['./externalModules/GradedVoronoiDiagram/AdaptiveGraphGenerator.py ./out/StressField_Tet_v2.stress ', ...
				sprintf('%g ',latticeSizeCtrl), sprintf('%g',sizeAspectRatio)];
			pyrunfile(callGradedVoronoiGenerater_python);
			LoadGeneratedGraphFromFileObj_B();
			
			%Voxelization%
			if optEdgeAlignmentComparison_
				[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveEles] = MGD_VoxelizeMeshEdges_PerEdge(edgeWidth, passiveElements);
			else
				MGD_VoxelizeMeshEdges();
				voxelsAlongLatticeEdges = MGD_CoatMeshEdgesWithVoxels_B(edgeWidth, passiveElements);				
			end
			volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
			disp(['............Determining Lower Bound of Lattice Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
				sprintf(' with Size Ctrl Para %g', latticeSizeCtrl)]);	
			if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
				densityLayout_(voxelsAlongLatticeEdges) = 1; return;
			end				
			if volumeFractionDesign_ < targetDepositionRatio
				upperLatticeSizeCtrl = lowerLatticeSizeCtrl;
				lowerLatticeSizeCtrl = lowerLatticeSizeCtrl / 1.5;
			else
				break;
			end		
			if lowerLatticeSizeCtrl < 0.001
				warning('Inappropriate settings for the material budget!');
				densityLayout_(voxelsAlongLatticeEdges) = 1;
				return;
			end			
		end	
	end

	%%Determine the target Lattice Size density control
	idx = 1;
	while abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio > permittedVolumeDeviation			
		latticeSizeCtrl = (lowerLatticeSizeCtrl + upperLatticeSizeCtrl) / 2;

		%%Create Graded Voronoi Diagram
		callGradedVoronoiGenerater_python = ['./externalModules/GradedVoronoiDiagram/AdaptiveGraphGenerator.py ./out/StressField_Tet_v2.stress ', ...
			sprintf('%g ',latticeSizeCtrl), sprintf('%g',sizeAspectRatio)];
		pyrunfile(callGradedVoronoiGenerater_python);
		LoadGeneratedGraphFromFileObj_B();
		
		if optEdgeAlignmentComparison_
			[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveEles] = MGD_VoxelizeMeshEdges_PerEdge(edgeWidth, passiveElements);
		else
			MGD_VoxelizeMeshEdges();
			voxelsAlongLatticeEdges = MGD_CoatMeshEdgesWithVoxels_B(edgeWidth, passiveElements);				
		end
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Design Iteration ', sprintf('%d', idx), sprintf('. Design Volume Fraction: %.6f', volumeFractionDesign_), ...
			sprintf(' with Line Density Para %g', latticeSizeCtrl)]);
		
		if volumeFractionDesign_>targetDepositionRatio		
			lowerLatticeSizeCtrl = latticeSizeCtrl;
		else
			upperLatticeSizeCtrl = latticeSizeCtrl;
		end
		idx = idx + 1;
		if idx > 10
			warning('Stress-aware Graded Voronoi Diagram Infill failed to converge to the prescribed design'); break;
		end
	end	
	densityLayout_(voxelsAlongLatticeEdges) = 1;
	voxelsOnBoundary_ = setdiff(voxelsOnBoundary_, voxelsAlongLatticeEdgesWithoutPassiveEles);
	voxelsOnBoundary_ = setdiff(voxelsOnBoundary_, unique([voxelsInLoadingArea_(:); voxelsInFixingArea_(:)]));
	tEnd = toc(tStart);
	disp(['............Conduct Stress-aligned Conforming Lattice Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);	
end

function LoadGeneratedGraphFromFileObj_B()
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	
	%%Read Field-aligned Graph in
	IO_ImportVertexEdgeGraph('./out/StressField_Tet_v2_Voronoi.obj');	
	frameStruct4Voxelization_ = vertexEdgeGraph_;
	frameStruct4Voxelization_.edgeLengths = vecnorm(frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,1),:) ...
		- frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,2),:),2,2);
end

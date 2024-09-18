function SAGS_StressAwareGradedVoronoiDiagramGeneration(edgeWidth, targetDepositionRatio, numLayerboundary, ...
                numLayerLoads, numLayerFixation, sizeAspectRatio)
	global meshHierarchy_;
	global volumeFractionDesign_; 
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;	
	global densityLayout_;
	global densityLayout4Vis_;
	
	upperLatticeSizeCtrl = 0.1;
	lowerLatticeSizeCtrl = 0.04;
	
	permittedVolumeDeviation = 0.05;
	opt_DetermingLowerBound = 1;
	densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	densityLayout4Vis_ = zeros(size(meshHierarchy_(1).eleMapForward));
	if targetDepositionRatio>0.9
		warning('Close to a solid domain, no need for design!');
		densityLayout_ = ones(size(densityLayout_));
		densityLayout4Vis_(meshHierarchy_(1).eleMapBack) = 1;
		volumeFractionDesign_ = 1;
		tEnd = toc(tStart);
		disp(['............Conduct Stress-aware Graded Voronoi Diagram Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);			
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
		densityLayout4Vis_(meshHierarchy_(1).eleMapBack(passiveElements),1) = 1;
		tEnd = toc(tStart);
		disp(['............Conduct Stress-aligned Conforming Lattice Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);		
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
		[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);
		
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Determining Upper Bound of Voronoi Cell Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Size Ctrl Para %g', latticeSizeCtrl)]);
		if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
			densityLayout_(voxelsAlongLatticeEdges) = 1; 
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
			densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;
			tEnd = toc(tStart);
			disp(['............Conduct Stress-aware Graded Voronoi Diagram Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);			
			return;
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
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
			densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;			
			tEnd = toc(tStart);
			disp(['............Conduct Stress-aware Graded Voronoi Diagram Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);			
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
			[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);
			
			volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
			disp(['............Determining Lower Bound of Voronoi Cell Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
				sprintf(' with Size Ctrl Para %g', latticeSizeCtrl)]);	
			if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
				densityLayout_(voxelsAlongLatticeEdges) = 1; 
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
				densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;				
				tEnd = toc(tStart);
				disp(['............Conduct Stress-aware Graded Voronoi Diagram Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);
				return;
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
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
				densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;				
				tEnd = toc(tStart);
				disp(['............Conduct Stress-aware Graded Voronoi Diagram Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);					
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
		
		[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);
		
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Design Iteration ', sprintf('%d', idx), sprintf('. Design Volume Fraction: %.6f', volumeFractionDesign_), ...
			sprintf(' with Voronoi Cell Size Para %g', latticeSizeCtrl)]);
		
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
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
	densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;	
	tEnd = toc(tStart);
	disp(['............Conduct Stress-aware Graded Voronoi Diagram Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);	
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

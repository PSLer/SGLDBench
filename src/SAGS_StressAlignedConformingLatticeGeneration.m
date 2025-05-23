function SAGS_StressAlignedConformingLatticeGeneration(edgeWidth, targetDepositionRatio, numLayerboundary, ...
                numLayerLoads, numLayerFixation, aspectRatio)
	global outPath_;
	global meshHierarchy_;
	global volumeFractionDesign_; 
	global voxelsOnBoundary_;
	global voxelsInLoadingArea_;
	global voxelsInFixingArea_;	
	global densityLayout_;
	global densityLayout4Vis_;
	global dataPrep4SAGS_;
    
    upperLatticeSizeCtrl = 0.8;  
    lowerLatticeSizeCtrl = 0.4;
	  
	permittedVolumeDeviation = 0.03;
	opt_DetermingLowerBound = 1;
	densityLayout_ = zeros(meshHierarchy_(1).numElements,1);
	densityLayout4Vis_ = zeros(size(meshHierarchy_(1).eleMapForward));
	if targetDepositionRatio>0.9
		warning('Close to a solid domain, no need for design!');
		densityLayout_ = ones(size(densityLayout_));
		densityLayout4Vis_(meshHierarchy_(1).eleMapBack) = 1;
		volumeFractionDesign_ = 1;
		tEnd = toc(tStart);
		disp(['............Conduct Stress-aligned Conforming Lattice Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);			
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
	
	%% Compute Frame Field (Principal Stress Field) from Stress Field
	
	%% Prepare Tet-Mesh for running "tensor-field-meshing.exe"
	fid = fopen(strcat(outPath_, 'FrameData4Gao2017.mesh'), 'w');	
	fprintf(fid, '%d %d\n', [size(dataPrep4SAGS_.nodeCoords,1) size(dataPrep4SAGS_.eNodMat,1)]);
	fprintf(fid, '%.6f %.6f %.6f\n', dataPrep4SAGS_.nodeCoords');
	fprintf(fid, '%d %d %d %d %d\n', [4*ones(size(dataPrep4SAGS_.eNodMat,1),1) dataPrep4SAGS_.eNodMat-1]');
	fclose(fid);	
	
    if exist('../externalModules/Gao2017/', 'dir')
	    callGao2017_Executable = strcat('"../externalModules/Gao2017/tensor-field-meshing.exe" -b -i', char(strcat(" ", strcat(outPath_, 'FrameData4Gao2017'))));
    else
        callGao2017_Executable = strcat('"./externalModules/Gao2017/tensor-field-meshing.exe" -b -i', char(strcat(" ", strcat(outPath_, 'FrameData4Gao2017'))));
    end
	%% Determine the upper bound for lattice size control
	volumeFractionDesign_ = 1;
	while volumeFractionDesign_ > targetDepositionRatio
		latticeSizeCtrl = upperLatticeSizeCtrl;
		SetupFrameFieldFile(latticeSizeCtrl, aspectRatio);

		%%Run Gao2017
		system(callGao2017_Executable); pause(1);
		LoadGeneratedGraphFromFileObj();
		
		%Voxelization%
		[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);
		
		volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
		disp(['............Determining Upper Bound of Lattice Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
			sprintf(' with Size Ctrl Para %g', latticeSizeCtrl)]);
		if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
			densityLayout_(voxelsAlongLatticeEdges) = 1; 
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
			densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;
			return;
		end				
		if volumeFractionDesign_ > targetDepositionRatio
			lowerLatticeSizeCtrl = upperLatticeSizeCtrl; opt_DetermingLowerBound = 0;
			upperLatticeSizeCtrl = upperLatticeSizeCtrl * 1.25;
		else
			break;
		end
		if upperLatticeSizeCtrl > 4
			warning('Inappropriate settings for the material budget!');
			densityLayout_(voxelsAlongLatticeEdges) = 1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
			densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
			densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;			
			tEnd = toc(tStart);
			disp(['............Conduct Stress-aligned Conforming Lattice Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);				
			return;
		end		
	end

	if 	opt_DetermingLowerBound
		%% Determine the lower bound for lattice size control
		volumeFractionDesign_ = 0;
		while volumeFractionDesign_ < targetDepositionRatio
			latticeSizeCtrl = lowerLatticeSizeCtrl;
			SetupFrameFieldFile(latticeSizeCtrl, aspectRatio);
	
			%%Run Gao2017
			system(callGao2017_Executable); pause(1);
			LoadGeneratedGraphFromFileObj();
			
			%Voxelization%
			[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);
			
			volumeFractionDesign_ = numel(voxelsAlongLatticeEdges) / meshHierarchy_(1).numElements;
			disp(['............Determining Lower Bound of Lattice Size Control: ', sprintf('Volume Fraction %.6f', volumeFractionDesign_), ...
				sprintf(' with Size Ctrl Para %g', latticeSizeCtrl)]);		
			if abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio <= permittedVolumeDeviation
				densityLayout_(voxelsAlongLatticeEdges) = 1;
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
				densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;					
				tEnd = toc(tStart);
				disp(['............Conduct Stress-aligned Conforming Lattice Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);					
				return;
			end	
			if volumeFractionDesign_ < targetDepositionRatio
				upperLatticeSizeCtrl = lowerLatticeSizeCtrl;
				lowerLatticeSizeCtrl = lowerLatticeSizeCtrl / 1.25;
			else
				break;
			end		
			if lowerLatticeSizeCtrl < 0.3
				warning('Inappropriate settings for the material budget! Please increase the resolution of gateway tet-mesh!');
				densityLayout_(voxelsAlongLatticeEdges) = 1;
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
				densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
				densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;			
				tEnd = toc(tStart);
				disp(['............Conduct Stress-aligned Conforming Lattice Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);					
				return;
			end			
		end	
	end

	%%Determine the target Lattice Size density control
	idx = 1;
	while abs(volumeFractionDesign_-targetDepositionRatio) / targetDepositionRatio > permittedVolumeDeviation			
		latticeSizeCtrl = (lowerLatticeSizeCtrl + upperLatticeSizeCtrl) / 2;
		SetupFrameFieldFile(latticeSizeCtrl, aspectRatio);

		%%Run Gao2017
		system(callGao2017_Executable); pause(1);
		LoadGeneratedGraphFromFileObj();
		
		[voxelsAlongLatticeEdges, voxelsAlongLatticeEdgesWithoutPassiveElesMapback] = MGD_VoxelizeMeshEdges_PerEdge_B(edgeWidth, passiveElements);
		
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
			warning('Stress-aligned Conforming Lattice Infill failed to converge to the prescribed design'); break;
		end
	end	
	densityLayout_(voxelsAlongLatticeEdges) = 1;
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack(voxelsOnBoundary_),1) = -1;
	densityLayout4Vis_(meshHierarchy_(1).eleMapBack([voxelsInFixingArea_(:); voxelsInLoadingArea_(:)]),1) = 1;
	densityLayout4Vis_(voxelsAlongLatticeEdgesWithoutPassiveElesMapback) = 1;
	tEnd = toc(tStart);
	disp(['............Conduct Stress-aligned Conforming Lattice Infill Design Costs: ', sprintf('%.1f', tEnd), 's']);	
end

function SetupFrameFieldFile(latticeSizeCtrl, aspectRatio)
	global dataPrep4SAGS_;
	global outPath_;

	aspectRatio = min([1.0, aspectRatio]);
if 0 %%anisotropic cell
	%% Initialize Frame Field	
	dataPrep4SAGS_.frameField = dataPrep4SAGS_.ps(:, [2 3 4 6 7 8 10 11 12 1 5 9]);
	if 1==aspectRatio
		dataPrep4SAGS_.frameField(:,[10 11 12]) = latticeSizeCtrl * ones(size(dataPrep4SAGS_.frameField,1),3);
	else
		minVal = aspectRatio;
		dataPrep4SAGS_.frameField(:,[10 11 12]) = abs(dataPrep4SAGS_.frameField(:,[10 11 12]));
		for ii=1:size(dataPrep4SAGS_.nodeCoords,1)
			iFrame = abs(dataPrep4SAGS_.frameField(ii,[10 11 12]));
			iFrame = iFrame / max(iFrame);
			iFrame = max(iFrame, minVal);
			dataPrep4SAGS_.frameField(ii,[10 11 12]) = latticeSizeCtrl * iFrame;
		end	
	end
else %%anisotropic cell but varying sizes
	dataPrep4SAGS_.frameField = dataPrep4SAGS_.ps(:, [2 3 4 6 7 8 10 11 12 1 5 9]);	
    vonMises = max(dataPrep4SAGS_.vM, 1.0e-12);
    vonMises = vonMises.^(1/2);
    upperBound = 1/aspectRatio; lowerBound = 1.0;
    scalingFactors = (vonMises - min(vonMises)) / (max(vonMises)-min(vonMises)) * (upperBound-lowerBound) + lowerBound;
    scalingFactors = 1 ./scalingFactors;
	dataPrep4SAGS_.frameField(:,[10 11 12]) = latticeSizeCtrl * ones(size(dataPrep4SAGS_.frameField,1),3) .* scalingFactors;
	
end	
	dataPrep4SAGS_.frameField = dataPrep4SAGS_.frameField';
	dataPrep4SAGS_.frameField = dataPrep4SAGS_.frameField(:);
	dataPrep4SAGS_.frameField = reshape(dataPrep4SAGS_.frameField, 3, 4*size(dataPrep4SAGS_.nodeCoords,1))';
	%%Write Frame Field
	fid = fopen(strcat(outPath_, 'FrameData4Gao2017.txt'), 'w');	
	fprintf(fid, '%.6f %.6f %.6f\n', dataPrep4SAGS_.frameField');
	fclose(fid);	
end

function LoadGeneratedGraphFromFileObj()
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	global outPath_;

	%%Read Field-aligned Graph in
	IO_ImportVertexEdgeGraph(strcat(outPath_, 'FrameData4Gao2017_graph_opt.obj'));	
	frameStruct4Voxelization_ = vertexEdgeGraph_;
	frameStruct4Voxelization_.edgeLengths = vecnorm(frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,1),:) ...
		- frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,2),:),2,2);
end
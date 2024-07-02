function SAGS_SubdivideHexIntoTets()
	global meshHierarchy_;
	global inputSolidMesh_;
	global cellSize_;
	global modulus_; 
	global poissonRatio_; 
	global complianceSolid_;
	
	solidMesh4Sim = inputSolidMesh_;
	
	%Align Solid Mesh 2 Cartesian Mesh, Transfering BCs
	refBoundingBox_ = [min(meshHierarchy_(1).boundaryNodeCoords,[],1); max(meshHierarchy_(1).boundaryNodeCoords,[],1)];
	newOrigin = refBoundingBox_(1,:);
	boundingBoxSimMesh = [min(solidMesh4Sim.nodeCoords,[],1); max(solidMesh4Sim.nodeCoords,[],1)];
	solidMesh4Sim.nodeCoords = solidMesh4Sim.nodeCoords + (newOrigin - boundingBoxSimMesh(1,:));
	newCharacterDimension = max(refBoundingBox_(2,:)-refBoundingBox_(1,:))*1; %%  	
	boundingBoxSimMesh = [min(solidMesh4Sim.nodeCoords, [], 1); max(solidMesh4Sim.nodeCoords, [], 1)];
	solidMesh4Sim.nodeCoords = boundingBoxSimMesh(1,:) + (solidMesh4Sim.nodeCoords - boundingBoxSimMesh(1,:)) ...
		* (newCharacterDimension/max(boundingBoxSimMesh(2,:)-boundingBoxSimMesh(1,:)));
	boundingBoxSimMesh = [min(solidMesh4Sim.nodeCoords, [], 1); max(solidMesh4Sim.nodeCoords, [], 1)];
	solidMesh4Sim.boundaryNodeCoords = solidMesh4Sim.nodeCoords(1==solidMesh4Sim.nodeState,:);
	[loads, fixations] = TransferBCsFromVoxels2ArbitraryMesh(solidMesh4Sim);
	newCharacterDimension = max(refBoundingBox_(2,:)-refBoundingBox_(1,:))*cellSize_;
	solidMesh4Sim.nodeCoords = boundingBoxSimMesh(1,:) + (solidMesh4Sim.nodeCoords - boundingBoxSimMesh(1,:)) ...
		* (newCharacterDimension/max(boundingBoxSimMesh(2,:)-boundingBoxSimMesh(1,:)));
	boundingBoxSimMesh = [min(solidMesh4Sim.nodeCoords, [], 1); max(solidMesh4Sim.nodeCoords, [], 1)];
	solidMesh4Sim.boundaryNodeCoords = solidMesh4Sim.nodeCoords(1==solidMesh4Sim.nodeState,:);	
	
	%FEA with pure HEX or TET mesh
	[deformation, cartesianStress, compliance] = FEA_SolidTetOrHexMesh(solidMesh4Sim, loads, fixations, modulus_, poissonRatio_);
	if 1
		disp(['Compliance (Voxel): ', sprintf('%f', complianceSolid_)]);
		disp(['Compliance (Hex/Tet): ', sprintf('%f', compliance)]);	
	end

	% WrapFEAmodel(solidMesh4Sim, loads, fixations);
	
	%%Subdivision
	dataPrep4Gao2017 = struct('nodeCoords', [], 'eNodMat', [], 'cartesianStress', [], 'frame', []);
	if strcmp(solidMesh4Sim.meshType, 'HEX')
		[hex_to_tet_nodeCoords, hex_to_tet_eNodMat, hex_to_tet_cartesianStress] = SAGS_ConvertHex2Tet12(solidMesh4Sim, cartesianStress);
		dataPrep4Gao2017.nodeCoords = hex_to_tet_nodeCoords;
		dataPrep4Gao2017.eNodMat = hex_to_tet_eNodMat;
		dataPrep4Gao2017.cartesianStress = hex_to_tet_cartesianStress;
	else
		dataPrep4Gao2017.nodeCoords = solidMesh4Sim.nodeCoords;
		dataPrep4Gao2017.eNodMat = solidMesh4Sim.eNodMat;
		dataPrep4Gao2017.cartesianStress = cartesianStress;		
	end
end

function [loads, fixations] = TransferBCsFromVoxels2ArbitraryMesh(solidMesh4Sim)
	global meshHierarchy_;
	global loadingCond_;
	global fixingCond_;
	
	boundaryNodesArbitraryMesh = find(1==solidMesh4Sim.nodeState);
	numBoundaryNodes = numel(boundaryNodesArbitraryMesh);
	%%Loading
	coordsBasedLoads = [meshHierarchy_(1).boundaryNodeCoords(loadingCond_(:,1),:) loadingCond_(:,2:end)];
	allBoundaryNodes = zeros(numBoundaryNodes,1);
	loads = zeros(numBoundaryNodes, 1+3); %%Solid
	for ii=1:size(coordsBasedLoads,1)
		iLoad = coordsBasedLoads(ii,:);
		iCoord = iLoad(1:3);
		iForce = iLoad(4:end);
		[~, tarNode] = min(vecnorm(iCoord-solidMesh4Sim.boundaryNodeCoords, 2, 2));
		loads(tarNode,1) = tarNode;
		loads(tarNode,2:end) = loads(tarNode,2:end) + iForce;
		allBoundaryNodes(tarNode) = allBoundaryNodes(tarNode) + 1;
	end
	realLoadedNodes = find(allBoundaryNodes);
	loads = loads(realLoadedNodes,:);
	loads(:,1) = boundaryNodesArbitraryMesh(loads(:,1));
	
	%%Fixation
	coordsBasedFixations = [meshHierarchy_(1).boundaryNodeCoords(fixingCond_(:,1),:) fixingCond_(:,2:end)];
	allBoundaryNodes = zeros(numBoundaryNodes,1);
	fixations = zeros(numBoundaryNodes, 1+3); %%Solid
	for ii=1:size(coordsBasedFixations,1)
		iFC = coordsBasedFixations(ii,:);
		iCoord = iFC(1:3);
		iFixation = iFC(4:end);
		[~, tarNode] = min(vecnorm(iCoord-solidMesh4Sim.boundaryNodeCoords, 2, 2));
		fixations(tarNode,1) = tarNode;
		fixations(tarNode,2:end) = iFixation;
		allBoundaryNodes(tarNode) = allBoundaryNodes(tarNode) + 1;
	end
	realFixedNodes = find(allBoundaryNodes);
	fixations = fixations(realFixedNodes,:);
	fixations(:,1) = boundaryNodesArbitraryMesh(fixations(:,1));
end

function WrapFEAmodel(solidMesh4Sim, loads, fixations)
	materialIndicatorField_ = ones(solidMesh4Sim.numElements,1);
	fileName = './out/Data4MiniFEM.MiniFEM';
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s ', 'Version');
	fprintf(fid, '%.1f\n', 2.0);
	switch solidMesh4Sim.meshType			
		case 'TET'
			fprintf(fid, '%s %s ', 'Solid Tet');
			fprintf(fid, '%d\n', 1);
			
			fprintf(fid, '%s ', 'Vertices:');
			fprintf(fid, '%d\n', solidMesh4Sim.numNodes);		
			fprintf(fid, '%.6e %.6e %.6e\n', solidMesh4Sim.nodeCoords');

			fprintf(fid, '%s ', 'Elements:');
			fprintf(fid, '%d \n', solidMesh4Sim.numElements);
			fprintf(fid, '%d %d %d %d %d\n', [solidMesh4Sim.eNodMat materialIndicatorField_]');

			fprintf(fid, '%s %s ', 'Node Forces:'); 
			fprintf(fid, '%d\n', size(loads,1));
			if ~isempty(loads)
				fprintf(fid, '%d %.6e %.6e %.6e\n', loads');
			end
			fprintf(fid, '%s %s ', 'Fixed Nodes:'); fprintf(fid, '%d\n', size(fixations,1));
			if ~isempty(fixations)
				fprintf(fid, '%d %d %d %d\n', fixations');						
			end		
		case 'HEX'
			fprintf(fid, '%s %s ', 'Solid Hex');
			fprintf(fid, '%d\n', 1);
			
			fprintf(fid, '%s ', 'Vertices:');
			fprintf(fid, '%d\n', solidMesh4Sim.numNodes);		
			fprintf(fid, '%.6e %.6e %.6e\n', solidMesh4Sim.nodeCoords');

			fprintf(fid, '%s ', 'Elements:');
			fprintf(fid, '%d \n', solidMesh4Sim.numElements);
			fprintf(fid, '%d %d %d %d %d %d %d %d %d\n', [solidMesh4Sim.eNodMat materialIndicatorField_]');

			fprintf(fid, '%s %s ', 'Node Forces:'); 
			fprintf(fid, '%d\n', size(loads,1));
			if ~isempty(loads)
				fprintf(fid, '%d %.6e %.6e %.6e\n', loads');
			end
			fprintf(fid, '%s %s ', 'Fixed Nodes:'); fprintf(fid, '%d\n', size(fixations,1));
			if ~isempty(fixations)
				fprintf(fid, '%d %d %d %d\n', fixations');						
			end
	end
	fclose(fid);	
end

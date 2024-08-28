function SAGS_ObtainStressField_ArbitrarySolidMesh()
	global meshHierarchy_;
	global inputSolidMesh_;
	global cellSize_;
	global modulus_; 
	global poissonRatio_; 
	global complianceSolid_;
	global dataPrep4SAGS_;
	
	solidMesh4Sim = inputSolidMesh_;
	
	%Align Solid Mesh 2 Cartesian Mesh, Transfering BCs
	refBoundingBox_ = [min(meshHierarchy_(1).boundaryNodeCoords,[],1); max(meshHierarchy_(1).boundaryNodeCoords,[],1)];
	newOrigin = refBoundingBox_(1,:);
	boundingBoxSimMesh = [min(solidMesh4Sim.nodeCoords,[],1); max(solidMesh4Sim.nodeCoords,[],1)];
	solidMesh4Sim.nodeCoords = solidMesh4Sim.nodeCoords + (newOrigin - boundingBoxSimMesh(1,:));
	newCharacterDimension = max(refBoundingBox_(2,:)-refBoundingBox_(1,:))*cellSize_; %%  	
	boundingBoxSimMesh = [min(solidMesh4Sim.nodeCoords, [], 1); max(solidMesh4Sim.nodeCoords, [], 1)];
	solidMesh4Sim.nodeCoords = boundingBoxSimMesh(1,:) + (solidMesh4Sim.nodeCoords - boundingBoxSimMesh(1,:)) ...
		* (newCharacterDimension/max(boundingBoxSimMesh(2,:)-boundingBoxSimMesh(1,:)));
	boundingBoxSimMesh = [min(solidMesh4Sim.nodeCoords, [], 1); max(solidMesh4Sim.nodeCoords, [], 1)];
	solidMesh4Sim.boundaryNodeCoords = solidMesh4Sim.nodeCoords(1==solidMesh4Sim.nodeState,:);
	[loads, fixations] = Common_TransferBCsFromVoxels2ArbitraryMesh(solidMesh4Sim);
	
	%FEA with pure HEX or TET mesh
	[deformation, cartesianStress, compliance] = FEA_SolidTetOrHexMesh(solidMesh4Sim, loads, fixations, modulus_, poissonRatio_);
	if 1
		disp(['Compliance (Voxel): ', sprintf('%f', complianceSolid_)]);
		disp(['Compliance (Hex/Tet): ', sprintf('%f', compliance)]);	
	end

	% WrapFEAmodel(solidMesh4Sim, loads, fixations);
	
	%%Subdivision
	dataPrep4SAGS_ = struct('nodeCoords', [], 'eNodMat', [], 'cartesianStress', [], 'ps', [], 'frameField', []);
	if strcmp(solidMesh4Sim.meshType, 'HEX')
		[hex_to_tet_nodeCoords, hex_to_tet_eNodMat, hex_to_tet_cartesianStress] = SAGS_ConvertHex2Tet12(solidMesh4Sim, cartesianStress);
		dataPrep4SAGS_.nodeCoords = hex_to_tet_nodeCoords;
		dataPrep4SAGS_.eNodMat = hex_to_tet_eNodMat;
		dataPrep4SAGS_.cartesianStress = hex_to_tet_cartesianStress;
	else
		dataPrep4SAGS_.nodeCoords = solidMesh4Sim.nodeCoords;
		dataPrep4SAGS_.eNodMat = solidMesh4Sim.eNodMat;
		dataPrep4SAGS_.cartesianStress = cartesianStress;		
	end
	dataPrep4SAGS_.ps = zeros(size(dataPrep4SAGS_.nodeCoords,1),12);
	for ii=1:size(dataPrep4SAGS_.nodeCoords,1)
		dataPrep4SAGS_.ps(ii,:)= FEA_ComputePrincipalStress(dataPrep4SAGS_.cartesianStress(ii,:));
	end
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

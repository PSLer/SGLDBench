function SAGS_AdaptCartesianStressData()
	global meshHierarchy_;
	global boundingBox_;
	global cartesianStressField_;
	global dataPrep4SAGS_;
	
	solidMesh4Sim = struct('nodeCoords', [], 'eNodMat', [], 'numElements', [], 'numNodes', []); 
	solidMesh4Sim.eNodMat = meshHierarchy_(1).eNodMat;
	solidMesh4Sim.numElements = meshHierarchy_(1).numElements;
	solidMesh4Sim.numNodes = meshHierarchy_(1).numNodes;
	
	nx = meshHierarchy_(1).resX;
    ny = meshHierarchy_(1).resY;
    nz = meshHierarchy_(1).resZ;
	xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nx:boundingBox_(2,1);
	ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/ny:boundingBox_(1,2);
	zSeed = boundingBox_(1,3):(boundingBox_(2,3)-boundingBox_(1,3))/nz:boundingBox_(2,3);
	solidMesh4Sim.nodeCoords = zeros(meshHierarchy_(1).numNodes,3);
	tmp = repmat(reshape(repmat(xSeed,ny+1,1), (nx+1)*(ny+1), 1), (nz+1), 1);
	nodPosX = reshape(tmp, ny+1, nx+1, nz+1);
	solidMesh4Sim.nodeCoords(:,1) = tmp(meshHierarchy_(1).nodMapBack,1);
	tmp = repmat(repmat(ySeed,1,nx+1)', (nz+1), 1);
	nodPosY = reshape(tmp, ny+1, nx+1, nz+1);
	solidMesh4Sim.nodeCoords(:,2) = tmp(meshHierarchy_(1).nodMapBack,1);
	tmp = reshape(repmat(zSeed,(nx+1)*(ny+1),1), (nx+1)*(ny+1)*(nz+1), 1);		
	nodPosZ = reshape(tmp, ny+1, nx+1, nz+1);
	solidMesh4Sim.nodeCoords(:,3) = tmp(meshHierarchy_(1).nodMapBack,1);
	
	[hex_to_tet_nodeCoords, hex_to_tet_eNodMat, hex_to_tet_cartesianStress] = SAGS_ConvertHex2Tet12(solidMesh4Sim, cartesianStressField_);
	dataPrep4SAGS_ = struct('nodeCoords', [], 'eNodMat', [], 'cartesianStress', [], 'ps', [], 'frameField', []);
	dataPrep4SAGS_.nodeCoords = hex_to_tet_nodeCoords;
	dataPrep4SAGS_.eNodMat = hex_to_tet_eNodMat;
	dataPrep4SAGS_.cartesianStress = hex_to_tet_cartesianStress;
	
	dataPrep4SAGS_.ps = zeros(size(dataPrep4SAGS_.nodeCoords,1),12);
	for ii=1:size(dataPrep4SAGS_.nodeCoords,1)
		dataPrep4SAGS_.ps(ii,:)= FEA_ComputePrincipalStress(dataPrep4SAGS_.cartesianStress(ii,:));
	end
end
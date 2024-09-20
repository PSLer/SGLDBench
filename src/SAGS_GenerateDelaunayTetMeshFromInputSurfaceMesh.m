function SAGS_GenerateDelaunayTetMeshFromInputSurfaceMesh(targetNumberOfTetElements)
	global outPath_;
	global surfaceTriMesh_;
	global gateWayTetMesh_;
	global meshHierarchy_;
	
	%%Align Dimensions to Voxel Model
	gateWayTetMesh_ = surfaceTriMesh_;
	refBoundingBox = [min(meshHierarchy_(1).boundaryNodeCoords,[],1); max(meshHierarchy_(1).boundaryNodeCoords,[],1)];
% refBoundingBox	
	newOrigin = refBoundingBox(1,:);
	newCharacterDimension = max(refBoundingBox(2,:)-refBoundingBox(1,:));
	boundingBoxGatewayMesh = [min(gateWayTetMesh_.nodeCoords,[],1); max(gateWayTetMesh_.nodeCoords,[],1)];
	gateWayTetMesh_.nodeCoords = gateWayTetMesh_.nodeCoords + (newOrigin - boundingBoxGatewayMesh(1,:));
	boundingBoxGatewayMesh = [min(gateWayTetMesh_.nodeCoords, [], 1); max(gateWayTetMesh_.nodeCoords, [], 1)];
	gateWayTetMesh_.nodeCoords = boundingBoxGatewayMesh(1,:) + (gateWayTetMesh_.nodeCoords - boundingBoxGatewayMesh(1,:)) ...
		* (newCharacterDimension/max(boundingBoxGatewayMesh(2,:)-boundingBoxGatewayMesh(1,:)));
	boundingBoxGatewayMesh = [min(gateWayTetMesh_.nodeCoords, [], 1); max(gateWayTetMesh_.nodeCoords, [], 1)];
% boundingBoxGatewayMesh	
    gateWayTetMesh_.nodeCoords(:,1) = gateWayTetMesh_.nodeCoords(:,1) + refBoundingBox(2,1)-boundingBoxGatewayMesh(2,1);
    gateWayTetMesh_.nodeCoords(:,2) = gateWayTetMesh_.nodeCoords(:,2) + refBoundingBox(2,2)-boundingBoxGatewayMesh(2,2);
    gateWayTetMesh_.nodeCoords(:,3) = gateWayTetMesh_.nodeCoords(:,3) + refBoundingBox(2,3)-boundingBoxGatewayMesh(2,3);
	
	%%Write out Gateway Triangular Surface Mesh for Tet-meshing by tetgen
	IO_ExportTriSurfMesh_PLY(gateWayTetMesh_, strcat(outPath_, 'gatewayMesh.ply'));
	
	%%Running TetGen
	totVolume = meshHierarchy_(1).eleSize(1) * meshHierarchy_(1).eleSize(2) * meshHierarchy_(1).eleSize(3) * meshHierarchy_(1).numElements;
	% targetNumberOfTetElements = 100000;
	maxCellVol = round(totVolume/targetNumberOfTetElements);
    if exist('../externalModules/TetGen/', 'dir')
        callTetGen = strcat('"../externalModules/TetGen/tetgen.exe" -gq1.414a', num2str(maxCellVol), char(strcat(" ", strcat(outPath_, 'gatewayMesh.ply'))));
    else
        callTetGen = strcat('"./externalModules/TetGen/tetgen.exe" -gq1.414a', num2str(maxCellVol), char(strcat(" ", strcat(outPath_, 'gatewayMesh.ply'))));
    end
	
% callTetGen
	disp('Generating Tet-mesh with tetgen ...');
	system(callTetGen);
	
	%%Load Gateway Tet-mesh
	gateWayTetMesh_ = IO_ImportSolidTetMesh_MESH_TetGen(strcat(outPath_, 'gatewayMesh.1.mesh'));	
end
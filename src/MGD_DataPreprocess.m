%%Mesh/Graph based Structural Design
function MGD_DataPreprocess()
    global inputSolidMesh_;
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	global meshHierarchy_;
	
	%%Align Dimensions
	refBoundingBox_ = [min(meshHierarchy_(1).boundaryNodeCoords,[],1); max(meshHierarchy_(1).boundaryNodeCoords,[],1)];  
	newOrigin = refBoundingBox_(1,:);
	newCharacterDimension = max(refBoundingBox_(2,:)-refBoundingBox_(1,:));
	boundingBoxFrame = [min(frameStruct4Voxelization_.nodeCoords,[],1); max(frameStruct4Voxelization_.nodeCoords,[],1)];
	frameStruct4Voxelization_.nodeCoords = frameStruct4Voxelization_.nodeCoords + (newOrigin - boundingBoxFrame(1,:));
	boundingBoxFrame = [min(frameStruct4Voxelization_.nodeCoords, [], 1); max(frameStruct4Voxelization_.nodeCoords, [], 1)];
	frameStruct4Voxelization_.nodeCoords = boundingBoxFrame(1,:) + (frameStruct4Voxelization_.nodeCoords - boundingBoxFrame(1,:)) ...
		* (newCharacterDimension/max(boundingBoxFrame(2,:)-boundingBoxFrame(1,:)));
	boundingBoxFrame = [min(frameStruct4Voxelization_.nodeCoords, [], 1); max(frameStruct4Voxelization_.nodeCoords, [], 1)];
    frameStruct4Voxelization_.nodeCoords(:,1) = frameStruct4Voxelization_.nodeCoords(:,1) + refBoundingBox_(2,1)-boundingBoxFrame(2,1);
    frameStruct4Voxelization_.nodeCoords(:,2) = frameStruct4Voxelization_.nodeCoords(:,2) + refBoundingBox_(2,2)-boundingBoxFrame(2,2);
    frameStruct4Voxelization_.nodeCoords(:,3) = frameStruct4Voxelization_.nodeCoords(:,3) + refBoundingBox_(2,3)-boundingBoxFrame(2,3);
	boundingBoxFrame = [min(frameStruct4Voxelization_.nodeCoords, [], 1); max(frameStruct4Voxelization_.nodeCoords, [], 1)]; 
	frameStruct4Voxelization_.edgeLengths = vecnorm(frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,1),:) ...
		- frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,2),:),2,2);	
end
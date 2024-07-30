%%Mesh/Graph based Structural Design
function MGD_DataPreprocess()
    global inputSolidMesh_;
	global vertexEdgeGraph_;
	global frameStruct4Voxelization_;
	global meshHierarchy_;
	
	%%Extract Graph for voxelization
	frameStruct4Voxelization_ = vertexEdgeGraph_;
	if ~vertexEdgeGraph_.state
		switch inputSolidMesh_.meshType
			case 'HEX'
				cell2edge = [1 2  2 3  3 4  4 1  5 6  6 7  7 8  8 5  1 5  2 6  3 7  4 8];
				numEdgePerCell = 12;
			case 'TET'
				cell2edge = [1 2  2 3  3 1  1 4  2 4  3 4];
				numEdgePerCell = 6;					
		end
		edgeIndices = inputSolidMesh_.eNodMat(:, cell2edge)';
		edgeIndices = reshape(edgeIndices, 2, numEdgePerCell*inputSolidMesh_.numElements);
		tmp = sort(edgeIndices,1)';
		[frameStruct4Voxelization_.eNodMat, ~, ~] = unique(tmp, 'rows');
		frameStruct4Voxelization_.numEdges = size(frameStruct4Voxelization_.eNodMat,1);
		frameStruct4Voxelization_.nodeCoords = inputSolidMesh_.nodeCoords;
		frameStruct4Voxelization_.numNodes = size(frameStruct4Voxelization_.nodeCoords,1);	
	end
	
	%%Align Dimensions
	refBoundingBox_ = [min(meshHierarchy_(1).boundaryNodeCoords,[],1); max(meshHierarchy_(1).boundaryNodeCoords,[],1)];
refBoundingBox_    
	newOrigin = refBoundingBox_(1,:);
	newCharacterDimension = max(refBoundingBox_(2,:)-refBoundingBox_(1,:));
	boundingBoxFrame = [min(frameStruct4Voxelization_.nodeCoords,[],1); max(frameStruct4Voxelization_.nodeCoords,[],1)];
	frameStruct4Voxelization_.nodeCoords = frameStruct4Voxelization_.nodeCoords + (newOrigin - boundingBoxFrame(1,:));
	boundingBoxFrame = [min(frameStruct4Voxelization_.nodeCoords, [], 1); max(frameStruct4Voxelization_.nodeCoords, [], 1)];
	frameStruct4Voxelization_.nodeCoords = boundingBoxFrame(1,:) + (frameStruct4Voxelization_.nodeCoords - boundingBoxFrame(1,:)) ...
		* (newCharacterDimension/max(boundingBoxFrame(2,:)-boundingBoxFrame(1,:)));
    frameStruct4Voxelization_.nodeCoords(:,1) = frameStruct4Voxelization_.nodeCoords(:,1) + refBoundingBox_(2,1)-boundingBoxFrame(2,1);
    frameStruct4Voxelization_.nodeCoords(:,2) = frameStruct4Voxelization_.nodeCoords(:,2) + refBoundingBox_(2,2)-boundingBoxFrame(2,2);
    frameStruct4Voxelization_.nodeCoords(:,3) = frameStruct4Voxelization_.nodeCoords(:,3) + refBoundingBox_(2,3)-boundingBoxFrame(2,3);
	boundingBoxFrame = [min(frameStruct4Voxelization_.nodeCoords, [], 1); max(frameStruct4Voxelization_.nodeCoords, [], 1)];
boundingBoxFrame    
	frameStruct4Voxelization_.edgeLengths = vecnorm(frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,1),:) ...
		- frameStruct4Voxelization_.nodeCoords(frameStruct4Voxelization_.eNodMat(:,2),:),2,2);	
end
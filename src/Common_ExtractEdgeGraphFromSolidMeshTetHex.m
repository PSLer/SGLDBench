function edgeGraph = Common_ExtractEdgeGraphFromSolidMeshTetHex(solidMesh)
	
	edgeGraph = Data_VertexEdgeGraphStruct();
	switch solidMesh.meshType
		case 'HEX'
			cell2edge = [1 2  2 3  3 4  4 1  5 6  6 7  7 8  8 5  1 5  2 6  3 7  4 8];
			numEdgePerCell = 12;
		case 'TET'
			cell2edge = [1 2  2 3  3 1  1 4  2 4  3 4];
			numEdgePerCell = 6;					
	end
	
	edgeIndices = solidMesh.eNodMat(:, cell2edge)';
	edgeIndices = reshape(edgeIndices, 2, numEdgePerCell*solidMesh.numElements);
	tmp = sort(edgeIndices,1)';
	[edgeGraph.eNodMat, ~, ~] = unique(tmp, 'rows');
	edgeGraph.numEdges = size(edgeGraph.eNodMat,1);
	edgeGraph.nodeCoords = solidMesh.nodeCoords;
	edgeGraph.numNodes = size(edgeGraph.nodeCoords,1);
end
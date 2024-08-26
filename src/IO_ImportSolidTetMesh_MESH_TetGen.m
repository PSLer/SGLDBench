function gatewayTetMesh = IO_ImportSolidTetMesh_MESH_TetGen(fileName)
	
	gatewayTetMesh = Data_ArbitraryMeshStruct();
	
	%%head
	fid = fopen(fileName, 'r');
	while 1
        tmp = fscanf(fid, '%s', 1);
        if strcmp(tmp, 'Vertices'), break; end
    end
    %%Nodes
	numNodes = fscanf(fid, '%d', 1);
	nodeCoords = fscanf(fid, '%f %f %f  %f', [4, numNodes]); 
	nodeCoords = nodeCoords'; nodeCoords(:,4) = [];
    
	%%Edges
    fileEdges = fscanf(fid, '%s', 1);
    numEdges = fscanf(fid, '%d', 1);
    edges = fscanf(fid, '%d %d %d', [3, numEdges]); 
    edges = edges'; edges(:,3) = [];
    
    %%Triangles
    while 1
        tmp = fscanf(fid, '%s', 1);
        if strcmp(tmp, 'Triangles'), break; end
    end
    tmp = fscanf(fid, '%s', 1);
    numTriangles = fscanf(fid, '%d', 1);
    triFaces = fscanf(fid, '%d %d %d %d', [4, numTriangles]); 
    
	%%Tets
    while 1
        tmp = fscanf(fid, '%s', 1);
        if strcmp(tmp, 'Tetrahedra'), break; end
    end
    tmp = fscanf(fid, '%s', 1);
	numEles = fscanf(fid, '%d', 1);
	eNodMat = fscanf(fid, '%d %d %d %d %d', [5, numEles]); 
	eNodMat = eNodMat'; eNodMat(:,end) = [];

	fclose(fid);
	
	gatewayTetMesh.state = 1;
	gatewayTetMesh.numNodes = numNodes;
	gatewayTetMesh.nodeCoords = nodeCoords;
	gatewayTetMesh.numElements = numEles;
	gatewayTetMesh.eNodMat = eNodMat;
end
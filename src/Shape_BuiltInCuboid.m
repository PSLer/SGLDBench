function Shape_BuiltInCuboid(xDim, yDim, zDim)
	global surfaceTriMesh_;
	
	surfaceTriMesh_.numNodes = 8;
	surfaceTriMesh_.nodeCoords = [
		0		0		0
		xDim	0		0
		xDim	yDim	0
		0		yDim	0
		0		0		zDim
		xDim	0		zDim
		xDim	yDim	zDim
		0		yDim	zDim		
	];
	surfaceTriMesh_.numElements = 12;
	surfaceTriMesh_.eNodMat = [
		1 	2	3
		3	4	1
		5	8	7
		7	6	5
		2	6	7
		7	3	2
		8	5	1
		1	4	8
		1	5	6
		2	1	6
		3	7	8
		8	4	3
	];
	surfaceTriMesh_.state = 1;
end
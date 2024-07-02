function Shape_BuiltInLshape(xDim, yDim, zDim, cPointX, cPointZ)
	global surfaceTriMesh_;
	
	surfaceTriMesh_.numNodes = 12;
	surfaceTriMesh_.nodeCoords = [
		0		0		0
		xDim	0		0
		xDim	yDim	0
		0		yDim	0
		0		0		zDim
		cPointX	0		zDim
		cPointX	yDim	zDim
		0		yDim	zDim
		cPointX	0		cPointZ
		xDim	0		cPointZ
		xDim	yDim	cPointZ
		cPointX	yDim	cPointZ
	
	];
	surfaceTriMesh_.numElements = 20;
	surfaceTriMesh_.eNodMat = [
		1 	2	3
		3	4	1
		5	8	7
		7	6	5
		9	12	11
		11	10	9
		8	5	1
		1	4	8
		6	7	12
		12	9	6
		2	10	11
		11	3	2
		9	5	6
		9	1	5
		9	2	1
		9	10	2
		12	7	8
		12	8	4
		12	4	3
		12	3	11
	];
	surfaceTriMesh_.state = 1;
end
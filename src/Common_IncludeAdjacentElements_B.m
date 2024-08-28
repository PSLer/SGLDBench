function oEleList = Common_IncludeAdjacentElements_B(iEleList)
	global meshHierarchy_;
	global eleX_;
	global eleY_; 
	global eleZ_;
	
	iEleListMapBack = meshHierarchy_(1).eleMapBack(iEleList);
	%%	1	4	7		10	 13	  16		19	 22	  25
	%%	2	5	8		11	 14*  17		20	 23	  26
	%%	3	6	9		12	 15   18		21	 24	  27
	%%	 bottom				middle				top
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;
	% [eleX, eleY, eleZ] = Common_NodalizeDesignDomain([resX-1 resY-1 resZ-1], [1 1 1; resX resY resZ]);
	eleX = eleX_(iEleListMapBack);
	eleY = eleY_(iEleListMapBack);
	eleZ = eleZ_(iEleListMapBack);
	
	tmpX = [eleX-1 eleX-1 eleX-1  eleX eleX eleX  eleX+1 eleX+1 eleX+1];
	tmpX = [tmpX tmpX tmpX]; tmpX = tmpX(:);
	tmpY = [eleY+1 eleY eleY-1  eleY+1 eleY eleY-1  eleY+1 eleY eleY-1]; 
	tmpY = [tmpY tmpY tmpY]; tmpY = tmpY(:);
	tmpZ = [eleZ eleZ eleZ eleZ eleZ eleZ eleZ eleZ eleZ];
	tmpZ = [tmpZ-1 tmpZ tmpZ+1]; tmpZ = tmpZ(:);
	xNegative = find(tmpX<1); xPositive = find(tmpX>resX);
	yNegative = find(tmpY<1); yPositive = find(tmpY>resY);
	zNegative = find(tmpZ<1); zPositive = find(tmpZ>resZ);
	allInvalidEles = unique([xNegative; xPositive; yNegative; yPositive; zNegative; zPositive]);
	tmpX(allInvalidEles) = []; tmpY(allInvalidEles) = []; tmpZ(allInvalidEles) = [];
	oEleListMapBack = resX*resY*(tmpZ-1) + resY*(tmpX-1) + resY-tmpY + 1;
	oEleList = meshHierarchy_(1).eleMapForward(oEleListMapBack);
	oEleList(oEleList<1) = []; oEleList = unique(oEleList);
end
function [voxelsAlongLatticeEdgesAndPassiveElements, voxelsAlongLatticeEdges] = MGD_VoxelizeMeshEdges_PerEdge_B(numLayersVoxelAroundEdge, passiveElements)
	global meshHierarchy_;
	global frameStruct4Voxelization_; 
	global eleX_;
	global eleY_; 
	global eleZ_;
	
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;
	[eleX_, eleY_, eleZ_] = Common_NodalizeDesignDomain([resX-1 resY-1 resZ-1], [1 1 1; resX resY resZ]);

	%% Edge Voxelization
	%%Distribute Sampling Points Along Mesh Edges
	refEleSize = meshHierarchy_(1).eleSize(1) * sqrt(2)/2; %% 8/9 or sqrt(2)/2 is just a random scaling factor
	maxEdgeLength = max(frameStruct4Voxelization_.edgeLengths);
	testNumEdges = frameStruct4Voxelization_.numEdges; %% frameStruct4Voxelization_.numEdges, 2000
	assoVoxels = struct('iEdge', 0, 'arr', []);
	assoVoxels = repmat(assoVoxels, testNumEdges, 1);
	% 
	for ii=1:testNumEdges
		iEdge = frameStruct4Voxelization_.eNodMat(ii,:);
		iEdgeEndPoints = frameStruct4Voxelization_.nodeCoords(iEdge,:);
		iNumSamplingPoints = ceil(frameStruct4Voxelization_.edgeLengths(ii) / refEleSize) + 1;
		iSamplingPoints = zeros(iNumSamplingPoints,3);
		iSamplingPoints(:,1) = linspace(iEdgeEndPoints(1,1), iEdgeEndPoints(2,1), iNumSamplingPoints);
		iSamplingPoints(:,2) = linspace(iEdgeEndPoints(1,2), iEdgeEndPoints(2,2), iNumSamplingPoints);
		iSamplingPoints(:,3) = linspace(iEdgeEndPoints(1,3), iEdgeEndPoints(2,3), iNumSamplingPoints);
		assoVoxels(ii).arr = zeros(1,iNumSamplingPoints);
		for jj=1:iNumSamplingPoints
			jPoint = iSamplingPoints(jj,:);	
			[eleIndex, opt] = LocatePointOnCartesianMesh_B(jPoint);
			if opt
				assoVoxels(ii).arr(jj) = eleIndex;
			end
		end
		assoVoxels(ii).arr(assoVoxels(ii).arr<1) = [];
		if isempty(assoVoxels(ii).arr), continue; end
		assoVoxels(ii).arr = unique(assoVoxels(ii).arr);
		for kk=1:numLayersVoxelAroundEdge-1
			assoVoxels(ii).arr = Common_IncludeAdjacentElements_C(assoVoxels(ii).arr);
		end
		assoVoxels(ii).arr = assoVoxels(ii).arr(:)';
		assoVoxels(ii).iEdge = ii;
	end
	emptyEdges = [];
	for ii=1:testNumEdges
		if isempty(assoVoxels(ii).arr)
			emptyEdges(end+1,1) = ii;
		end
	end
	assoVoxels(emptyEdges) = [];
	voxelsAlongLatticeEdges = unique([assoVoxels.arr])';
	for ii=1:numel(assoVoxels)
		iEles = assoVoxels(ii).arr;
		iEles = meshHierarchy_(1).eleMapForward(iEles);
		iEles(iEles<1) = [];
		assoVoxels(ii).arr = iEles(:)';
	end	
	voxelsAlongLatticeEdgesMapforward = unique([assoVoxels.arr])';
	voxelsAlongLatticeEdgesAndPassiveElements = unique([voxelsAlongLatticeEdgesMapforward(:); passiveElements(:)]);
	frameStruct4Voxelization_.AssociatedVoxels = assoVoxels;
end

function [eleIndex, opt] = LocatePointOnCartesianMesh_B(physicalCoordinates)
	global meshHierarchy_;
	global boundingBox_;
	
	eleIndex = 0; opt = 0;
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;
	physicalCoordinates = physicalCoordinates - boundingBox_(1,:);
	if 0==physicalCoordinates(1)
		eleX = 1;				
	else
		eleX = ceil(physicalCoordinates(1)/meshHierarchy_(1).eleSize(1));
		if eleX<1 || eleX>resX, return; end
	end
	if 0==physicalCoordinates(2)
		eleY = 1;
	else
		eleY = ceil(physicalCoordinates(2)/meshHierarchy_(1).eleSize(2));
		if eleY<1 || eleY>resY, return; end
	end
	if 0==physicalCoordinates(3)
		eleZ = 1;
	else
		eleZ = ceil(physicalCoordinates(3)/meshHierarchy_(1).eleSize(3));
		if eleZ<1 || eleZ>resZ, return; end
	end			
	
	eleIndex = resX*resY*(eleZ-1) + resY*(eleX-1)+(resY-eleY+1);
	%eleIndex = meshHierarchy_(1).eleMapForward(eleIndex);
	if eleIndex, opt = 1; end	
end

function oEleList = Common_IncludeAdjacentElements_C(iEleList)
	global meshHierarchy_;
	global eleX_;
	global eleY_; 
	global eleZ_;
	
	%iEleListMapBack = meshHierarchy_(1).eleMapBack(iEleList);
	iEleListMapBack = iEleList;
	%%	1	4	7		10	 13	  16		19	 22	  25
	%%	2	5	8		11	 14*  17		20	 23	  26
	%%	3	6	9		12	 15   18		21	 24	  27
	%%	 bottom				middle				top
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;
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
	% oEleList = meshHierarchy_(1).eleMapForward(oEleListMapBack);
	% oEleList(oEleList<1) = []; 
	oEleList = oEleListMapBack;
	oEleList = unique(oEleList);
end
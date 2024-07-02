function volumeFractionOfVoxelizedMeshEdges = MGD_VoxelizeMeshEdges()
	global meshHierarchy_;
	global frameStruct4Voxelization_; 
	global voxelizedMeshEdges_;
	global densityLayoutWithoutBoundary_;
	%% Edge Voxelization
	%%Distribute Sampling Points Along Mesh Edges
	refEleSize = meshHierarchy_(1).eleSize(1) * sqrt(2)/2; %% 8/9 or sqrt(2)/2 is just a random scaling factor
	maxEdgeLength = max(frameStruct4Voxelization_.edgeLengths);
	maxSamplingPointsPerEdge = ceil(maxEdgeLength/refEleSize)+1;
	samplingPointsOnMeshEdges = [];
	samplingPointsOnMeshEdges = zeros(maxSamplingPointsPerEdge*frameStruct4Voxelization_.numEdges,3);
	testNumEdges = frameStruct4Voxelization_.numEdges; %% frameStruct4Voxelization_.numEdges, 2000
	iSlider = 0;
	for ii=1:testNumEdges
		iEdge = frameStruct4Voxelization_.eNodMat(ii,:);
		iEdgeEndPoints = frameStruct4Voxelization_.nodeCoords(iEdge,:);
		iNumSamplingPoints = ceil(frameStruct4Voxelization_.edgeLengths(ii) / refEleSize) + 1;
		iSamplingPoints = zeros(iNumSamplingPoints,3);
		iSamplingPoints(:,1) = linspace(iEdgeEndPoints(1,1), iEdgeEndPoints(2,1), iNumSamplingPoints);
		iSamplingPoints(:,2) = linspace(iEdgeEndPoints(1,2), iEdgeEndPoints(2,2), iNumSamplingPoints);
		iSamplingPoints(:,3) = linspace(iEdgeEndPoints(1,3), iEdgeEndPoints(2,3), iNumSamplingPoints);
		% samplingPointsOnMeshEdges(end+1:end+iNumSamplingPoints,1:3) = iSamplingPoints;
		samplingPointsOnMeshEdges(iSlider+1:iSlider+iNumSamplingPoints,:) = iSamplingPoints;
		iSlider = iSlider + iNumSamplingPoints;
	end
	samplingPointsOnMeshEdges = samplingPointsOnMeshEdges(1:iSlider,:);

	%%Locate the Voxels at the Sampling Points

	numSamplingPointsAlongMeshEdges = size(samplingPointsOnMeshEdges,1);
	voxelizedMeshEdges_ = zeros(numSamplingPointsAlongMeshEdges,1);
	unrelatedSamplingPoints = [];
	for ii=1:numSamplingPointsAlongMeshEdges
		iSampingPoint = samplingPointsOnMeshEdges(ii,:);
		[eleIndex, opt] = LocatePointOnCartesianMesh(iSampingPoint);
		if opt
			voxelizedMeshEdges_(ii) = eleIndex;
		else
			% [~, eleIndexIdxOnBoundary] = min(vecnorm(iSampingPoint-boundaryElementsCentroidList,2,2));
			% voxelizedMeshEdges_(ii) = meshHierarchy_(1).elementsOnBoundary(eleIndexIdxOnBoundary);
			unrelatedSamplingPoints(end+1,1:3) = iSampingPoint;
		end
	end
	% aaaa = voxelizedMeshEdges_;
	voxelizedMeshEdges_(voxelizedMeshEdges_<1) = []; 
	% bbbb = voxelizedMeshEdges_;
	voxelizedMeshEdges_ = unique(voxelizedMeshEdges_);
	volumeFractionOfVoxelizedMeshEdges = numel(voxelizedMeshEdges_) / meshHierarchy_(1).numElements;
	densityLayoutWithoutBoundary_ = zeros(meshHierarchy_(1).numElements,1);
	densityLayoutWithoutBoundary_(voxelizedMeshEdges_(:)) = 1;
    totVolume = meshHierarchy_(1).numElements * (meshHierarchy_(1).eleSize(1)*meshHierarchy_(1).eleSize(2)*meshHierarchy_(1).eleSize(3));
	frameVolume = numel(voxelizedMeshEdges_) * (meshHierarchy_(1).eleSize(1)*meshHierarchy_(1).eleSize(2)*meshHierarchy_(1).eleSize(3));
end

function [eleIndex, opt] = LocatePointOnCartesianMesh(physicalCoordinates)
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
	
	tarEle = resX*resY*(eleZ-1) + resY*(eleX-1)+(resY-eleY+1);
	eleIndex = meshHierarchy_(1).eleMapForward(tarEle);
	if eleIndex, opt = 1; end	
end
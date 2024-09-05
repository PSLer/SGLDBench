function FEA_VoxelBasedDiscretization()	
	global outPath_;
	global nelx_; global nely_; global nelz_;
	global boundingBox_;
	global voxelizedVolume_; 
	global meshHierarchy_;
	global coarsestResolutionControl_;
	global eNodMatHalfTemp_;
	global numLevels_;
	%    z
	%    |__ x
	%   / 
	%  -y                            
	%            8--------------7      	
	%			/ |			   /|	
	%          5-------------6	|
	%          |  |          |  |
	%          |  |          |  |	
	%          |  |          |  |   
	%          |  4----------|--3  
	%     	   | /           | /
	%          1-------------2             
	%			Hexahedral element
	
	%%1. adjust voxel resolution for building Mesh Hierarchy
	[nely_, nelx_, nelz_] = size(voxelizedVolume_);
    numVoxels = numel(find(voxelizedVolume_));	
	if numVoxels<coarsestResolutionControl_
		error('There is no Sufficient Resolution to Build Mesh Hierarchy! Try Reducing the Coarsest Resolution or Increacing Original Resolution');
	end
	if isempty(numLevels_)
		numLevels_ = 0;
		while numVoxels>=coarsestResolutionControl_
			numLevels_ = numLevels_+1;
			numVoxels = round(numVoxels/8);
		end
	else
		numLevels_ = round(numLevels_)-1;
	end
	adjustedNelx = ceil(nelx_/2^numLevels_)*2^numLevels_;
	adjustedNely = ceil(nely_/2^numLevels_)*2^numLevels_;
	adjustedNelz = ceil(nelz_/2^numLevels_)*2^numLevels_;
	numLevels_ = numLevels_ + 1;
	if adjustedNelx>nelx_
		voxelizedVolume_(:,end+1:adjustedNelx,:) = false(nely_,adjustedNelx-nelx_,nelz_);
	end
	if adjustedNely>nely_
		voxelizedVolume_(end+1:adjustedNely,:,:) = false(adjustedNely-nely_,adjustedNelx,nelz_);
	end
	if adjustedNelz>nelz_
		voxelizedVolume_(:,:,end+1:adjustedNelz) = false(adjustedNely,adjustedNelx,adjustedNelz-nelz_);
	end

	%%2. initialize characteristic size
	boundingBox_ = [0 0 0; adjustedNelx adjustedNely adjustedNelz];
	
	%%3. initialize the finest mesh
	if 1==meshHierarchy_(1).state, meshHierarchy_ = Data_CartesianMeshStruct(); end
	meshHierarchy_.resX = adjustedNelx; nx = meshHierarchy_.resX;
	meshHierarchy_.resY = adjustedNely; ny = meshHierarchy_.resY;
	meshHierarchy_.resZ = adjustedNelz; nz = meshHierarchy_.resZ;
	meshHierarchy_.eleSize = (boundingBox_(2,:) - boundingBox_(1,:)) ./ [nx ny nz];

	%%4. identify solid&void elements
	voxelizedVolume_ = voxelizedVolume_(:);
	meshHierarchy_.eleMapBack = find(voxelizedVolume_);
	meshHierarchy_.eleMapBack = int32(meshHierarchy_.eleMapBack);
	meshHierarchy_.numElements = numel(meshHierarchy_.eleMapBack);
	meshHierarchy_.eleMapForward = zeros(nx*ny*nz,1, 'int32');	
	meshHierarchy_.eleMapForward(meshHierarchy_.eleMapBack) = (1:meshHierarchy_.numElements)';
	% meshHierarchy_.eleMapForward = int32(meshHierarchy_.eleMapForward);
		
	%%5. discretize
	nodenrs = reshape(1:(nx+1)*(ny+1)*(nz+1), 1+ny, 1+nx, 1+nz); nodenrs = int32(nodenrs);
	eNodVec = reshape(nodenrs(1:end-1,1:end-1,1:end-1)+1, nx*ny*nz, 1);
	eNodMatHalfTemp_ = repmat(eNodVec,1,8);
	eNodMat = repmat(eNodVec(meshHierarchy_.eleMapBack),1,8);	
	tmp = [0 ny+[1 0] -1 (ny+1)*(nx+1)+[0 ny+[1 0] -1]]; tmp = int32(tmp);
	for ii=1:8
		eNodMat(:,ii) = eNodMat(:,ii) + repmat(tmp(ii), meshHierarchy_.numElements,1);
		eNodMatHalfTemp_(:,ii) = eNodMatHalfTemp_(:,ii) + repmat(tmp(ii), nx*ny*nz,1);	
	end
	eNodMatHalfTemp_ = eNodMatHalfTemp_(:,[3 4 7 8]);
	meshHierarchy_.nodMapBack = unique(eNodMat);
	meshHierarchy_.numNodes = length(meshHierarchy_.nodMapBack);
	meshHierarchy_.numDOFs = meshHierarchy_.numNodes*3;
	meshHierarchy_.nodMapForward = zeros((nx+1)*(ny+1)*(nz+1),1, 'int32');
	meshHierarchy_.nodMapForward(meshHierarchy_.nodMapBack) = (1:meshHierarchy_.numNodes)';

	for ii=1:8
		eNodMat(:,ii) = meshHierarchy_.nodMapForward(eNodMat(:,ii));
	end	
	
	%%6. identify boundary info.	
	meshHierarchy_.numNod2ElesVec = zeros(meshHierarchy_.numNodes,1);
	for jj=1:8
		iNodes = eNodMat(:,jj);
		meshHierarchy_.numNod2ElesVec(iNodes) = meshHierarchy_.numNod2ElesVec(iNodes) + 1;
	end	

	meshHierarchy_.nodesOnBoundary = find(meshHierarchy_.numNod2ElesVec<8);
	meshHierarchy_.nodesOnBoundary = int32(meshHierarchy_.nodesOnBoundary);
	allNodes = zeros(meshHierarchy_.numNodes,1,'int32');
	allNodes(meshHierarchy_.nodesOnBoundary) = 1;	
	tmp = zeros(meshHierarchy_.numElements,1,'int32');
	for ii=1:8
		tmp = tmp + allNodes(eNodMat(:,ii));
	end
	meshHierarchy_.elementsOnBoundary = int32(find(tmp>0));
	blockIndex = Solving_MissionPartition(meshHierarchy_.numElements, 5.0e6);
	for ii=1:size(blockIndex,1)				
		rangeIndex = (blockIndex(ii,1):blockIndex(ii,2))';
		patchIndices = eNodMat(rangeIndex, [4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
		patchIndices = reshape(patchIndices(:), 4, 6*numel(rangeIndex));
		tmp = zeros(meshHierarchy_.numNodes, 1);
		tmp(meshHierarchy_.nodesOnBoundary) = 1;
		tmp = tmp(patchIndices); tmp = sum(tmp,1);
		iBoundaryEleFaces = patchIndices(:,find(4==tmp));
		meshHierarchy_.boundaryEleFaces(end+1:end+size(iBoundaryEleFaces,2),:) = iBoundaryEleFaces';
	end
	allNodes(meshHierarchy_.nodesOnBoundary) = (1:numel(meshHierarchy_.nodesOnBoundary))';
	meshHierarchy_.boundaryEleFaces = allNodes(meshHierarchy_.boundaryEleFaces);
	
	%%7. 
	% nodeCoords_ = zeros((nx+1)*(ny+1)*(nz+1),3);
	xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nx:boundingBox_(2,1); xSeed = single(xSeed);
	ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/ny:boundingBox_(1,2); ySeed = single(ySeed);
	zSeed = boundingBox_(1,3):(boundingBox_(2,3)-boundingBox_(1,3))/nz:boundingBox_(2,3); zSeed = single(zSeed);	
	eleCentroidList = zeros(meshHierarchy_.numElements,3, 'single');
	meshHierarchy_.boundaryNodeCoords = zeros(numel(meshHierarchy_(1).nodesOnBoundary),3, 'single');
	tmp = repmat(reshape(repmat(xSeed,ny+1,1), (nx+1)*(ny+1), 1), (nz+1), 1); tmp = tmp(meshHierarchy_.nodMapBack,1); %%Node Coords in X
	niftiwrite(tmp, strcat(outPath_, 'cache_coordX.nii'));
	meshHierarchy_.boundaryNodeCoords(:,1) = tmp(meshHierarchy_(1).nodesOnBoundary,1);
	for ii=1:size(blockIndex,1)	
		iSelEleNodes = eNodMat(blockIndex(ii,1):blockIndex(ii,2),:);
		eleCentX = tmp(iSelEleNodes);
		eleCentroidList(blockIndex(ii,1):blockIndex(ii,2),1) = sum(eleCentX,2)/8;
	end
	tmp = repmat(repmat(ySeed,1,nx+1)', (nz+1), 1); tmp = tmp(meshHierarchy_.nodMapBack,1); %%Node Coords in Y
	niftiwrite(tmp, strcat(outPath_, 'cache_coordY.nii'));
	meshHierarchy_.boundaryNodeCoords(:,2) = tmp(meshHierarchy_(1).nodesOnBoundary,1);
	for ii=1:size(blockIndex,1)	
		iSelEleNodes = eNodMat(blockIndex(ii,1):blockIndex(ii,2),:);
		eleCentY = tmp(iSelEleNodes);
		eleCentroidList(blockIndex(ii,1):blockIndex(ii,2),2) = sum(eleCentY,2)/8;
	end	
	tmp = reshape(repmat(zSeed,(nx+1)*(ny+1),1), (nx+1)*(ny+1)*(nz+1), 1); tmp = tmp(meshHierarchy_.nodMapBack,1); %Node Coords in Z
	niftiwrite(tmp, strcat(outPath_, 'cache_coordZ.nii'));
	meshHierarchy_.boundaryNodeCoords(:,3) = tmp(meshHierarchy_(1).nodesOnBoundary,1);
	for ii=1:size(blockIndex,1)	
		iSelEleNodes = eNodMat(blockIndex(ii,1):blockIndex(ii,2),:);
		eleCentZ = tmp(iSelEleNodes);
		eleCentroidList(blockIndex(ii,1):blockIndex(ii,2),3) = sum(eleCentZ,2)/8;
	end
	% meshHierarchy_.eNodMatHalf = eNodMat(:,[3 4 7 8]);
	meshHierarchy_.eNodMat = eNodMat;
	meshHierarchy_.state = 1;
	
	niftiwrite(eleCentroidList, strcat(outPath_, 'cache_eleCentroidList.nii'));
	% boundingBox_ = [min(meshHierarchy_.boundaryNodeCoords,[],1); max(meshHierarchy_.boundaryNodeCoords,[],1)];
end

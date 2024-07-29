function IO_ImportTopVoxels(fileName)
	global voxelizedVolume_;
	global nelx_; 
	global nely_; 
	global nelz_; 
	global fixingCond_; 
	global loadingCond_;
	global meshHierarchy_;
	global passiveElements_;	
	
	fid = fopen(fileName, 'r');
	fileHead = fscanf(fid, '%s %s %s %s', 4);
	if ~strcmp(fscanf(fid, '%s', 1), 'Resolution:'), error('Incompatible Mesh Data Format!'); end
	resolutions = fscanf(fid, '%d %d %d', 3);
	nelx_ = resolutions(1); nely_ = resolutions(2); nelz_ = resolutions(3);
	startReadSolidVoxels = fscanf(fid, '%s %s', 2);
	numSolidVoxels = fscanf(fid, '%d', 1);
	solidVoxels = fscanf(fid, '%d', [1 numSolidVoxels])';
	
	if ~strcmp(fscanf(fid, '%s', 1), 'Passive'), error('Incompatible Mesh Data Format!'); end
	if ~strcmp(fscanf(fid, '%s', 1), 'elements:'), error('Incompatible Mesh Data Format!'); end
	numPassiveElements = fscanf(fid, '%d', 1);
	passiveElements_ = fscanf(fid, '%d', [1 numPassiveElements])';
	if ~strcmp(fscanf(fid, '%s', 1), 'Fixations:'), error('Incompatible Mesh Data Format!'); end
	numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0		
		fixingCond_ = fscanf(fid, '%d %d %d %d', [4 numFixedNodes])';		
	end
	if ~strcmp(fscanf(fid, '%s', 1), 'Loads:'), error('Incompatible Mesh Data Format!'); end
	numLoadedNodes = fscanf(fid, '%d', 1);	
	if numLoadedNodes>0
		loadingCond_ = fscanf(fid, '%d %e %e %e', [4 numLoadedNodes])';								
	end
	fclose(fid);
	
	voxelizedVolume_ = false(nelx_*nely_*nelz_,1); 
	voxelizedVolume_(solidVoxels) = true;
	voxelizedVolume_ = reshape(voxelizedVolume_, nely_, nelx_, nelz_);
	FEA_VoxelBasedDiscretization();
	
	%%In case the resolution is slightly inconsistent	
	if ~isempty(fixingCond_)
		[~,uniqueFixedNodes] = unique(fixingCond_(:,1));
		fixingCond_ = fixingCond_(uniqueFixedNodes,:);
		[~,sortMap] = sort(fixingCond_(:,1),'ascend');
		fixingCond_ = fixingCond_(sortMap,:);
		fixingCond_ = AdaptBCExternalMdl(fixingCond_, [meshHierarchy_(1).resX+1 meshHierarchy_(1).resY+1 meshHierarchy_(1).resZ+1]);
		allNodes = zeros(size(meshHierarchy_(1).nodMapForward));
		allNodes(meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary)) = (1:numel(meshHierarchy_(1).nodesOnBoundary))';	
		fixingCond_(:,1) = allNodes(fixingCond_(:,1));		
	end
	if ~isempty(loadingCond_)
		[~,uniqueFixedNodes] = unique(loadingCond_(:,1));
		loadingCond_ = loadingCond_(uniqueFixedNodes,:);
		[~,sortMap] = sort(loadingCond_(:,1),'ascend');
		loadingCond_ = loadingCond_(sortMap,:);
		loadingCond_ = AdaptBCExternalMdl(loadingCond_, [meshHierarchy_(1).resX+1 meshHierarchy_(1).resY+1 meshHierarchy_(1).resZ+1]);
		allNodes = zeros(size(meshHierarchy_(1).nodMapForward));
		allNodes(meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary)) = (1:numel(meshHierarchy_(1).nodesOnBoundary))';	
		loadingCond_(:,1) = allNodes(loadingCond_(:,1));		
	end
	if ~isempty(passiveElements_)
		passiveElements_ = sort(passiveElements_, 'ascend');
		passiveElements_ = AdaptPassiveElementsExternalMdl(passiveElements_, [meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]);
		passiveElements_ = meshHierarchy_(1).eleMapForward(passiveElements_);
	end
end

function tarBC = AdaptBCExternalMdl(srcBC, adjustedRes)
	global nelx_; global nely_; global nelz_;
	nullNodeVolume = zeros((nelx_+1)*(nely_+1)*(nelz_+1),1);
	
	adjustedNnlx = adjustedRes(1);
	adjustedNnly = adjustedRes(2); 
	adjustedNnlz = adjustedRes(3);
	nnlx = nelx_+1;
	nnly = nely_+1;
	nnlz = nelz_+1;
	
	tmpBC = zeros(size(srcBC));
	[~, ia] = sort(srcBC(:,1));
	srcBC = srcBC(ia,:);
	nodeVolume = nullNodeVolume; 
	nodeVolume(srcBC(:,1)) = 1; 
	nodeVolume = reshape(nodeVolume, nely_+1, nelx_+1, nelz_+1);
	nodeVolume(:,end+1:adjustedNnlx,:) = zeros(nnly,adjustedNnlx-nnlx,nnlz);
	nodeVolume(end+1:adjustedNnly,:,:) = zeros(adjustedNnly-nnly,adjustedNnlx,nnlz);
	nodeVolume(:,:,end+1:adjustedNnlz) = zeros(adjustedNnly,adjustedNnlx,adjustedNnlz-nnlz);
	newLoadedNodes = find(nodeVolume);
	tmpBC(:,1) = newLoadedNodes;
	
	nodeVolume = nullNodeVolume; 
	nodeVolume(srcBC(:,1)) = srcBC(:,2); 
	nodeVolume = reshape(nodeVolume, nely_+1, nelx_+1, nelz_+1);
	nodeVolume(:,end+1:adjustedNnlx,:) = zeros(nnly,adjustedNnlx-nnlx,nnlz);
	nodeVolume(end+1:adjustedNnly,:,:) = zeros(adjustedNnly-nnly,adjustedNnlx,nnlz);
	nodeVolume(:,:,end+1:adjustedNnlz) = zeros(adjustedNnly,adjustedNnlx,adjustedNnlz-nnlz);
	nodeVolume = nodeVolume(:);
	tmpBC(:,2) = nodeVolume(newLoadedNodes);
	
	nodeVolume = nullNodeVolume; 
	nodeVolume(srcBC(:,1)) = srcBC(:,3); 
	nodeVolume = reshape(nodeVolume, nely_+1, nelx_+1, nelz_+1);
	nodeVolume(:,end+1:adjustedNnlx,:) = zeros(nnly,adjustedNnlx-nnlx,nnlz);
	nodeVolume(end+1:adjustedNnly,:,:) = zeros(adjustedNnly-nnly,adjustedNnlx,nnlz);
	nodeVolume(:,:,end+1:adjustedNnlz) = zeros(adjustedNnly,adjustedNnlx,adjustedNnlz-nnlz);
	nodeVolume = nodeVolume(:);
	tmpBC(:,3) = nodeVolume(newLoadedNodes);		
	
	nodeVolume = nullNodeVolume; 
	nodeVolume(srcBC(:,1)) = srcBC(:,4); 
	nodeVolume = reshape(nodeVolume, nely_+1, nelx_+1, nelz_+1);
	nodeVolume(:,end+1:adjustedNnlx,:) = zeros(nnly,adjustedNnlx-nnlx,nnlz);
	nodeVolume(end+1:adjustedNnly,:,:) = zeros(adjustedNnly-nnly,adjustedNnlx,nnlz);
	nodeVolume(:,:,end+1:adjustedNnlz) = zeros(adjustedNnly,adjustedNnlx,adjustedNnlz-nnlz);
	nodeVolume = nodeVolume(:);
	tmpBC(:,4) = nodeVolume(newLoadedNodes);			
	tarBC = tmpBC;
end

function adjustedVoxelIndices = AdaptPassiveElementsExternalMdl(srcElesMapback, adjustedRes)
	global nelx_; global nely_; global nelz_;
	nullVoxelVolume = zeros(nelx_*nely_*nelz_,1);
	
	adjustedNelx = adjustedRes(1);
	adjustedNely = adjustedRes(2); 
	adjustedNelz = adjustedRes(3);

	tmpBC = zeros(size(srcElesMapback));
	[~, ia] = sort(srcElesMapback(:,1));
	srcElesMapback = srcElesMapback(ia,:);
	voxelVolume = nullVoxelVolume; 
	voxelVolume(srcElesMapback(:,1)) = 1; 
	voxelVolume = reshape(voxelVolume, nely_, nelx_, nelz_);
	voxelVolume(:,end+1:adjustedNelx,:) = zeros(nely_,adjustedNelx-nelx_,nelz_);
	voxelVolume(end+1:adjustedNely,:,:) = zeros(adjustedNely-nely_,adjustedNelx,nelz_);
	voxelVolume(:,:,end+1:adjustedNelz) = zeros(adjustedNely,adjustedNelx,adjustedNelz-nelz_);
	adjustedVoxelIndices = find(voxelVolume);
end
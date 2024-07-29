function Temp_CreateFromWrappedVoxelFEAModel(fileName, varargin)
	%% CreateFromWrappedVoxelFEAmodel(fileName); 
	%% CreateFromWrappedVoxelFEAmodel(fileName, sizeScaling); 
	% global domainType_;
	global nelx_; 
	global nely_; 
	global nelz_; 
	global boundingBox_;
	global fixingCond_; 
	global loadingCond_;
	global meshHierarchy_;
    global voxelizedVolume_;
    
	fid = fopen(fileName, 'r');
	tmp = fscanf(fid, '%s %s', 2);
	domainType_ = fscanf(fid, '%s', 1);

	tmp = fscanf(fid, '%s', 1);
	tmp = fscanf(fid, '%d %d %d', 3);
	nelx_ = tmp(1); nely_ = tmp(2); nelz_ = tmp(3);
	if 1==nargin
		boundingBox_ = [0 0 0; nelx_ nely_ nelz_];
	else
		boundingBox_ = [0 0 0; varargin{1}*[nelx_ nely_ nelz_]/max([nelx_ nely_ nelz_])];
	end			
	tmp = fscanf(fid, '%s %s', 2);
	numEles = fscanf(fid, '%d', 1);
	eleList = fscanf(fid, '%d', [1 numEles])';
	allEles = false(nelx_*nely_*nelz_,1); allEles(eleList) = true;				
	tmp = fscanf(fid, '%s %s', 2);
	numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0		
		fixingCond_ = fscanf(fid, '%d %d %d %d', [4 numFixedNodes])';		
	end
	tmp = fscanf(fid, '%s %s', 2);
	numLoadedNodes = fscanf(fid, '%d', 1);	
	if numLoadedNodes>0
		loadingCond_ = fscanf(fid, '%d %e %e %e', [4 numLoadedNodes])';								
	end
	fclose(fid);
	
	voxelizedVolume_ = reshape(allEles, nely_, nelx_, nelz_);
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
function FEA_BuiltInBoundaryConditions4CuboidDesignDomain(opt)
	global boundingBox_;
	global nelx_; global nely_; global nelz_;
	global meshHierarchy_;
	global fixingCond_;
	global loadingCond_;
	DL = min(meshHierarchy_(1).boundaryNodeCoords,[],1);
	eleSize = meshHierarchy_(1).eleSize;
	switch opt
		case 'Cantilever 3'
			% 2.1.1 fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% 2.1.2 loading	
			coord = [nelx_ (2*meshHierarchy_(1).resY-nely_)/2 round(nelz_/2)] .* eleSize; force = [0.0 0.0 -1];
			[~,loadedNodes] = min(vecnorm(coord-meshHierarchy_(1).boundaryNodeCoords,2,2));
			loadingCond_ = [double(loadedNodes) force];
    		fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
		case 'Cantilever 4'
			% 2.1.1 fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% 2.1.2 loading
			coord = [nelx_ (2*meshHierarchy_(1).resY-nely_)/2 DL(3)] .* eleSize; force = [0.0 0.0 -1];
			[~,loadedNodes] = min(vecnorm(coord-meshHierarchy_(1).boundaryNodeCoords,2,2));
			loadingCond_ = [double(loadedNodes) force];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
		case 'Cantilever 5'
			% 2.1.1 fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% 2.1.2 loading
			coord = [nelx_ meshHierarchy_(1).resY-nely_ DL(3)] .* eleSize; force = [0.0 0.0 -1];
			[~,loadedNodes] = min(vecnorm(coord-meshHierarchy_(1).boundaryNodeCoords,2,2));
			loadingCond_ = [double(loadedNodes) force];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
		case 'Cantilever 6'
			% 2.1.1 fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% 2.1.2 loading
			loadedNodesTemp = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			loadedNodes = find(DL(3)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,3));
			loadedNodes = loadedNodesTemp(loadedNodes);
			numLoadedNodes = size(loadedNodes,1);
			force = [0.0 0.0 -1]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [double(loadedNodes) force];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
		case 'Cantilever 7'
			% 2.1.1 fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% 2.1.2 loading
			loadedNodes = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			numLoadedNodes = size(loadedNodes,1);
			force = [0.0 0.0 -1]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [double(loadedNodes) force];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
		case 'Cantilever 8'
			% 2.1.1 fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% 2.1.2 loading
			loadedNodesTemp = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			loadedNodes = find((meshHierarchy_(1).resY-nely_)*eleSize(2)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,2));
			loadedNodes = loadedNodesTemp(loadedNodes);
			numLoadedNodes = size(loadedNodes,1);
			force = [0.0 0.0 -1]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [double(loadedNodes) force];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)]; 
		case 'Cantilever 9'
			% 2.1.1 fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% 2.1.2 loading
			loadedNodesTemp = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			loadedNodes = find(round(nelz_/2)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,3));
			loadedNodes = loadedNodesTemp(loadedNodes);
			numLoadedNodes = size(loadedNodes,1);
			force = [1.0 0.0 0.0]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [double(loadedNodes) force];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
			
			loadedNodesTemp = find(DL(3)==meshHierarchy_(1).boundaryNodeCoords(:,3));
			loadedNodes = find(round(nelx_/2)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,1));
			loadedNodes = loadedNodesTemp(loadedNodes);
			
			force = [0.0 0.0 -1.0]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [loadingCond_; [double(loadedNodes) force]];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];			
	end
	
end
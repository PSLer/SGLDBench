function FEA_BuiltInBoundaryConditions4CuboidDesignDomain(opt)
	global boundingBox_;
	global nelx_; global nely_; global nelz_;
	global meshHierarchy_;
	global fixingCond_;
	global loadingCond_;
	DL = min(meshHierarchy_(1).boundaryNodeCoords,[],1);
	eleSize = meshHierarchy_(1).eleSize;
	switch opt
		case 'Cuboid - Cantilever 1'
			% fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
			% loading	
			coord = [nelx_ (2*meshHierarchy_(1).resY-nely_)/2 round(nelz_/2)] .* eleSize; force = [0.0 0.0 -1];
			[~,loadedNodes] = min(vecnorm(coord-meshHierarchy_(1).boundaryNodeCoords,2,2));
			loadingCond_ = [double(loadedNodes) force];   		
		case 'Cantilever 4'
			% fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			% loading
			coord = [nelx_ (2*meshHierarchy_(1).resY-nely_)/2 DL(3)] .* eleSize; force = [0.0 0.0 -1];
			[~,loadedNodes] = min(vecnorm(coord-meshHierarchy_(1).boundaryNodeCoords,2,2));
			loadingCond_ = [double(loadedNodes) force];
            fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
		case 'Cuboid - Cantilever 2'
			% fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));	
			fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
			% loading
			coord = [nelx_ meshHierarchy_(1).resY-nely_ DL(3)] .* eleSize; force = [0.0 0.0 -1];
			[~,loadedNodes] = min(vecnorm(coord-meshHierarchy_(1).boundaryNodeCoords,2,2));
			loadingCond_ = [double(loadedNodes) force];           
		case 'Cuboid - Cantilever 3'
			% fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
			% loading
			loadedNodesTemp = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			loadedNodes = find(DL(3)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,3));
			loadedNodes = loadedNodesTemp(loadedNodes);
			numLoadedNodes = size(loadedNodes,1);
			force = [0.0 0.0 -1]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [double(loadedNodes) force];          
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
		case 'Cuboid - Cantilever 4'
			% fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];			
			% loading
			loadedNodesTemp = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			loadedNodes = loadedNodesTemp(find(round(nelz_/2)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,3)));
			numLoadedNodes = size(loadedNodes,1);
			force = [1.0 0.0 0.0]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [double(loadedNodes) force];  
			
			loadedNodesTemp = find(DL(3)==meshHierarchy_(1).boundaryNodeCoords(:,3));
			loadedNodes = loadedNodesTemp(find(round(nelx_/2)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,1)));
			numLoadedNodes = size(loadedNodes,1);
			force = [0.0 0.0 -1.0]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [loadingCond_; [double(loadedNodes) force]];
		case 'Cuboid - MBB Half'
			% fixing
			fixedNodes = find(DL(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			fixingCond_ = [fixedNodes ones(length(fixedNodes),2) zeros(length(fixedNodes),1)]; 
			fixedNodesTemp = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			fixedNodes = fixedNodesTemp(find(DL(3)==meshHierarchy_(1).boundaryNodeCoords(fixedNodesTemp,3)));
			fixingCond_ = [fixingCond_; [fixedNodes ones(length(fixedNodes), 3)]];
			% loading
			loadedNodesTemp = find(0==meshHierarchy_(1).boundaryNodeCoords(:,1));
			loadedNodes = loadedNodesTemp(find(nelz_==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,3)));
			numLoadedNodes = size(loadedNodes,1);
			force = [0.0 0.0 -1]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [double(loadedNodes) force]; 
		case 'Cuboid - MBB'
			% fixing
			fixedNodesTemp = find(nelx_*eleSize(1)==meshHierarchy_(1).boundaryNodeCoords(:,1));
			fixedNodes = fixedNodesTemp(find(DL(3)==meshHierarchy_(1).boundaryNodeCoords(fixedNodesTemp,3)));
			fixingCond_ = [fixedNodes ones(length(fixedNodes), 3)];
			fixedNodesTemp = find(0==meshHierarchy_(1).boundaryNodeCoords(:,1));
			fixedNodes = fixedNodesTemp(find(DL(3)==meshHierarchy_(1).boundaryNodeCoords(fixedNodesTemp,3)));
			fixingCond_ = [fixingCond_; [fixedNodes ones(length(fixedNodes), 3)]];
			% loading
			loadedNodesTemp = find(nelz_*eleSize(3)==meshHierarchy_(1).boundaryNodeCoords(:,3));
			loadedNodes = loadedNodesTemp(find(round(nelx_/2)==meshHierarchy_(1).boundaryNodeCoords(loadedNodesTemp,1)));
			numLoadedNodes = size(loadedNodes,1);
			force = [0.0 0.0 -1.0]; iForce = force/numLoadedNodes; force = repmat(iForce,numLoadedNodes,1);
			loadingCond_ = [loadingCond_; [double(loadedNodes) force]];			
	end
	
end
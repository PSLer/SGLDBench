function [loads, fixations] = Common_TransferBCsFromVoxels2ArbitraryMesh(solidMesh4Sim)
	global meshHierarchy_;
	global loadingCond_;
	global fixingCond_;
	
	boundaryNodesArbitraryMesh = find(1==solidMesh4Sim.nodeState);
	numBoundaryNodes = numel(boundaryNodesArbitraryMesh);
	%%Loading
	coordsBasedLoads = [meshHierarchy_(1).boundaryNodeCoords(loadingCond_(:,1),:) loadingCond_(:,2:end)];
	allBoundaryNodes = zeros(numBoundaryNodes,1);
	loads = zeros(numBoundaryNodes, 1+3); %%Solid
	for ii=1:size(coordsBasedLoads,1)
		iLoad = coordsBasedLoads(ii,:);
		iCoord = iLoad(1:3);
		iForce = iLoad(4:end);
		[~, tarNode] = min(vecnorm(iCoord-solidMesh4Sim.boundaryNodeCoords, 2, 2));
		loads(tarNode,1) = tarNode;
		loads(tarNode,2:end) = loads(tarNode,2:end) + iForce;
		allBoundaryNodes(tarNode) = allBoundaryNodes(tarNode) + 1;
	end
	realLoadedNodes = find(allBoundaryNodes);
	loads = loads(realLoadedNodes,:);
	loads(:,1) = boundaryNodesArbitraryMesh(loads(:,1));
	
	%%Fixation
	coordsBasedFixations = [meshHierarchy_(1).boundaryNodeCoords(fixingCond_(:,1),:) fixingCond_(:,2:end)];
	allBoundaryNodes = zeros(numBoundaryNodes,1);
	fixations = zeros(numBoundaryNodes, 1+3); %%Solid
	for ii=1:size(coordsBasedFixations,1)
		iFC = coordsBasedFixations(ii,:);
		iCoord = iFC(1:3);
		iFixation = iFC(4:end);
		[~, tarNode] = min(vecnorm(iCoord-solidMesh4Sim.boundaryNodeCoords, 2, 2));
		fixations(tarNode,1) = tarNode;
		fixations(tarNode,2:end) = iFixation;
		allBoundaryNodes(tarNode) = allBoundaryNodes(tarNode) + 1;
	end
	realFixedNodes = find(allBoundaryNodes);
	fixations = fixations(realFixedNodes,:);
	fixations(:,1) = boundaryNodesArbitraryMesh(fixations(:,1));
end
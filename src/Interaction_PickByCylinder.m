function Interaction_PickByCylinder(axHandle, dir, varargin)
	global domainType_;
	global boundingBox_;
	global meshHierarchy_;
	global boundaryNodeCoords;
	global pickedNodeCache_;
	global hdPickedNode_;
	if strcmp(domainType_, '2D'), return; end
	
	
	if numel(pickedNodeCache_)<4
		warning('No Sufficient Nodes are Picked to Locate the Cylinder!'); return; 
	else
		ctr = sum(meshHierarchy_.boundaryNodeCoords(pickedNodeCache_(end-3:end),:),1) / 4;
		radi = sum(vecnorm(ctr-meshHierarchy_.boundaryNodeCoords(pickedNodeCache_(end-3:end),:),2,2)) / 4;
	end
	
	nodesWithinCircle = [];
	switch dir
		case 'X'
			if 2==nargin
				ctrlPoints = boundingBox_(:,1);
			else
				ctrlPoints = varargin{1};
			end
			ctrlPoints = sort(ctrlPoints, 'ascend');
			startP = ctrlPoints(1); endP = ctrlPoints(2);
			numLayers = (endP-startP) / meshHierarchy_(1).eleSize(1) + 1;	numLayers = round(numLayers);
			for ii=1:numLayers
				iLayerCtrlCoord = startP + (ii-1)*meshHierarchy_(1).eleSize(1);
				iTarBoundaryNodes = find(iLayerCtrlCoord == meshHierarchy_.boundaryNodeCoords(:,1));
				iSpottedBoundaryNodes = iTarBoundaryNodes(find(vecnorm(ctr([2 3]) - meshHierarchy_.boundaryNodeCoords(iTarBoundaryNodes,[2 3]),2,2)<=radi));
				if ~isempty(iSpottedBoundaryNodes)
					nodesWithinCircle(end+1:end+numel(iSpottedBoundaryNodes),1) = iSpottedBoundaryNodes;
				end
			end
		case 'Y'
			if 2==nargin
				ctrlPoints = boundingBox_(:,2);
			else
				ctrlPoints = varargin{1};
			end
			ctrlPoints = sort(ctrlPoints, 'ascend');			
			startP = ctrlPoints(1); endP = ctrlPoints(2);
			numLayers = (endP-startP) / meshHierarchy_(1).eleSize(2) + 1;	numLayers = round(numLayers);
			for ii=1:numLayers
				iLayerCtrlCoord = startP + (ii-1)*meshHierarchy_(1).eleSize(2);
				iTarBoundaryNodes = find(iLayerCtrlCoord == meshHierarchy_.boundaryNodeCoords(:,2));
				iSpottedBoundaryNodes = iTarBoundaryNodes(find(vecnorm(ctr([1 3]) - meshHierarchy_.boundaryNodeCoords(iTarBoundaryNodes,[1 3]),2,2)<=radi));
				if ~isempty(iSpottedBoundaryNodes)
					nodesWithinCircle(end+1:end+numel(iSpottedBoundaryNodes),1) = iSpottedBoundaryNodes;
				end
			end			
		case 'Z'
			if 2==nargin
				ctrlPoints = boundingBox_(:,3);
			else
				ctrlPoints = varargin{1};
			end
			ctrlPoints = sort(ctrlPoints, 'ascend');
			startP = ctrlPoints(1); endP = ctrlPoints(2);
			numLayers = (endP-startP) / meshHierarchy_(1).eleSize(3) + 1;	numLayers = round(numLayers);
			numLayers = (endP-startP) / meshHierarchy_(1).eleSize(2) + 1;	numLayers = round(numLayers);
			for ii=1:numLayers
				iLayerCtrlCoord = startP + (ii-1)*meshHierarchy_(1).eleSize(3);
				iTarBoundaryNodes = find(iLayerCtrlCoord == meshHierarchy_.boundaryNodeCoords(:,3));
				iSpottedBoundaryNodes = iTarBoundaryNodes(find(vecnorm(ctr([1 2]) - meshHierarchy_.boundaryNodeCoords(iTarBoundaryNodes,[1 2]),2,2)<=radi));
				if ~isempty(iSpottedBoundaryNodes)
					nodesWithinCircle(end+1:end+numel(iSpottedBoundaryNodes),1) = iSpottedBoundaryNodes;
				end
			end			
		otherwise
			warning('Wrong Input!'); return;
	end
	hold(axHandle, 'on');
	hdPickedNode_(end+1) = plot3(axHandle, meshHierarchy_.boundaryNodeCoords(nodesWithinCircle,1), meshHierarchy_.boundaryNodeCoords(nodesWithinCircle,2), ...
		meshHierarchy_.boundaryNodeCoords(nodesWithinCircle,3), '.r', 'LineWidth', 2, 'MarkerSize', 10);		
	numNewlyPickedNodes = length(nodesWithinCircle);
	pickedNodeCache_(end+1:end+numNewlyPickedNodes,1) = nodesWithinCircle;
end
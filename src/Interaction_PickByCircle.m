function Interaction_PickByCircle(axHandle, varargin)
	global meshHierarchy_;
	global boundaryNodeCoords;
	global pickedNodeCache_;
	global hdPickedNode_;
	
	if 1==nargin
		if size(pickedNodeCache_, 1) < 2
			warning('There MUST be at lease TWO Picked Nodes Available!'); return;
		end
		ctr = meshHierarchy_.boundaryNodeCoords(pickedNodeCache_(end-1),:);
		radi = norm(ctr-meshHierarchy_.boundaryNodeCoords(pickedNodeCache_(end),:));
		nodesWithinCircle = find(vecnorm(ctr-meshHierarchy_.boundaryNodeCoords,2,2)<=radi);
		nodesWithinCircle = setdiff(nodesWithinCircle, pickedNodeCache_(end-1:end));		
	elseif 2==nargin
		if isempty(pickedNodeCache_)
			warning('There MUST be at lease ONE Picked Nodes Available!'); return;
		end
		ctr = meshHierarchy_.boundaryNodeCoords(pickedNodeCache_(end),:);
		radi = varargin{1};
		nodesWithinCircle = find(vecnorm(ctr-meshHierarchy_.boundaryNodeCoords,2,2)<=radi);
		nodesWithinCircle = setdiff(nodesWithinCircle, pickedNodeCache_(end));			
	else
		warning('Wrong Input!'); return;
	end	
	hold(axHandle, 'on');
	hdPickedNode_(end+1) = plot3(axHandle, meshHierarchy_.boundaryNodeCoords(nodesWithinCircle,1), meshHierarchy_.boundaryNodeCoords(nodesWithinCircle,2), ...
		meshHierarchy_.boundaryNodeCoords(nodesWithinCircle,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);	
	numNewlyPickedNodes = length(nodesWithinCircle);
	pickedNodeCache_(end+1:end+numNewlyPickedNodes,1) = nodesWithinCircle;
end
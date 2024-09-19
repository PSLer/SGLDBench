function Interaction_PickBySelectionShpere(axHandle, sphereCtr, sphereRad)
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdPickedNode_;
	
	nodesWithinSelectionSphere = find(vecnorm(sphereCtr-meshHierarchy_(1).boundaryNodeCoords,2,2)<=sphereRad);
	
	% nodesWithinSelectionSphere = find(cP1(1)<=meshHierarchy_(1).boundaryNodeCoords(:,1));
	% nodesWithinSelectionSphere = nodesWithinSelectionSphere(find(cP2(1)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,1)));
	% nodesWithinSelectionSphere = nodesWithinSelectionSphere(find(cP1(2)<=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,2)));
	% nodesWithinSelectionSphere = nodesWithinSelectionSphere(find(cP2(2)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,2)));
	% nodesWithinSelectionSphere = nodesWithinSelectionSphere(find(cP1(3)<=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,3)));
	% nodesWithinSelectionSphere = nodesWithinSelectionSphere(find(cP2(3)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,3)));
	% nodesWithinSelectionSphere = setdiff(nodesWithinSelectionSphere, pickedNodeCache_);
	
	numNewlyPickedNodes = numel(nodesWithinSelectionSphere);
numNewlyPickedNodes	
	if numNewlyPickedNodes>0
		hold(axHandle, 'on');
		if isempty(hdPickedNode_)
			hdPickedNode_ = plot3(axHandle, meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,1), ...
				meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,2), ...
					meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);
		else
			hdPickedNode_(end+1) = plot3(axHandle, meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,1), ...
				meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,2), ...
					meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionSphere,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);		
		end
		pickedNodeCache_(end+1:end+numNewlyPickedNodes,1) = nodesWithinSelectionSphere;
	end
end
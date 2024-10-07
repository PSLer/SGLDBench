function Interaction_PickBySelectionBox(axHandle, cP1, cP2)
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdPickedNode_;

	for ii=1:3
		if cP1(ii)>cP2(ii), tmp = cP1(ii); cP1(ii) = cP2(ii); cP2(ii) = tmp; end
    end
	nodesWithinSelectionBox = find(cP1(1)<=meshHierarchy_(1).boundaryNodeCoords(:,1));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP2(1)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,1)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP1(2)<=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,2)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP2(2)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,2)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP1(3)<=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,3)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP2(3)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,3)));
	nodesWithinSelectionBox = setdiff(nodesWithinSelectionBox, pickedNodeCache_);
	
	numNewlyPickedNodes = numel(nodesWithinSelectionBox);
	if numNewlyPickedNodes>0
		hold(axHandle, 'on');
		if isempty(hdPickedNode_)
			hdPickedNode_ = plot3(axHandle, meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,1), ...
				meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,2), ...
					meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);
		else
			hdPickedNode_(end+1) = plot3(axHandle, meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,1), ...
				meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,2), ...
					meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);		
		end
		pickedNodeCache_(end+1:end+numNewlyPickedNodes,1) = nodesWithinSelectionBox;
	end
end
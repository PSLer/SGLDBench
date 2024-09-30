function Interaction_UnPickBySelectionBox(axHandle, cP1, cP2)
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdPickedNode_;
	
	if isempty(pickedNodeCache_), return; end
	
	nodesWithinSelectionBox = pickedNodeCache_(find(cP1(1)<=meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,1)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP2(1)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,1)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP1(2)<=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,2)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP2(2)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,2)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP1(3)<=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,3)));
	nodesWithinSelectionBox = nodesWithinSelectionBox(find(cP2(3)>=meshHierarchy_(1).boundaryNodeCoords(nodesWithinSelectionBox,3)));
	if isempty(nodesWithinSelectionBox), return; end
    % nodesWithinSelectionBox = setdiff(nodesWithinSelectionBox, pickedNodeCache_);
	
	set(hdPickedNode_, 'visible', 'off');
	pickedNodeCache_ = setdiff(pickedNodeCache_, nodesWithinSelectionBox);
	hold(axHandle, 'on');
    if isempty(pickedNodeCache_), return; end
	hdPickedNode_(end+1) = plot3(axHandle, meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,1), ...
		meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,2), ...
			meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);
end
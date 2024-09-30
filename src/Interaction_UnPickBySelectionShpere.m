function Interaction_UnPickBySelectionShpere(axHandle, sphereCtr, sphereRad)
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdPickedNode_;
	
	if isempty(pickedNodeCache_), return; end
	nodesWithinSelectionSphere = find(vecnorm(sphereCtr-meshHierarchy_(1).boundaryNodeCoords,2,2)<=sphereRad);

	if isempty(nodesWithinSelectionSphere), return; end
	
	set(hdPickedNode_, 'visible', 'off');
	pickedNodeCache_ = setdiff(pickedNodeCache_, nodesWithinSelectionSphere);
	hold(axHandle, 'on');
	hdPickedNode_(end+1) = plot3(axHandle, meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,1), ...
		meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,2), ...
			meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);
end
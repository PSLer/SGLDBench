function Interaction_PickByFace(axHandle, constantDir, opt)
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdPickedNode_;
	dcm_obj = datacursormode;
	info_struct = getCursorInfo(dcm_obj);
	if isempty(info_struct)
		warning('No Cursor Mode Available!'); return;
	end
	if ~(1==opt || 0==opt || -1==opt), warning('Wrong Input!'); return; end
	tarNode = info_struct.Position;
	hold(axHandle, 'on'); 
	switch constantDir
		case 'X'
			refNodePool = meshHierarchy_.boundaryNodeCoords(:,1);
			refPos = tarNode(1);				
		case 'Y'
			refNodePool = meshHierarchy_.boundaryNodeCoords(:,2);
			refPos = tarNode(2);				
		case 'Z'
			refNodePool = meshHierarchy_.boundaryNodeCoords(:,3);
			refPos = tarNode(3);				
		otherwise
			error('Wrongly Defined Line Direction!');
	end
	switch opt
		case 1
			nodesOnLine = find(refNodePool>=refPos);
		case 0
			nodesOnLine = find(abs(refNodePool-refPos)<=0.0);
		case -1
			nodesOnLine = find(refNodePool<=refPos);
	end	
	hdPickedNode_(end+1) = plot3(axHandle, meshHierarchy_.boundaryNodeCoords(nodesOnLine,1), meshHierarchy_.boundaryNodeCoords(nodesOnLine,2), ...
		meshHierarchy_.boundaryNodeCoords(nodesOnLine,3), 'xr', 'LineWidth', 2, 'MarkerSize', 10);
	numNewlyPickedNodes = length(nodesOnLine);
	pickedNodeCache_(end+1:end+numNewlyPickedNodes,1) = nodesOnLine;
end
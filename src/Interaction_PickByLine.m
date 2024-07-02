function Interaction_PickByLine(axHandle, lineDir, varargin)
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdpickedNode_;
	
	lineRelaxFactor = 0;
	dcm_obj = datacursormode;
	info_struct = getCursorInfo(dcm_obj);
	if isempty(info_struct)
		warning('No Cursor Mode Available!'); return;
	end
	tarNode = info_struct.Position;
	hold(axHandle, 'on');
	constantDir = ["X", "Y", "Z"];
	constantDir = setdiff(constantDir, lineDir);
	switch constantDir(1)
		case 'X'
			nodesOnLine1 = find(abs(meshHierarchy_.boundaryNodeCoords(:,1)-tarNode(1))<=lineRelaxFactor);
		case 'Y'
			nodesOnLine1 = find(abs(meshHierarchy_.boundaryNodeCoords(:,2)-tarNode(2))<=lineRelaxFactor);
		case 'Z'
			nodesOnLine1 = find(abs(meshHierarchy_.boundaryNodeCoords(:,3)-tarNode(3))<=lineRelaxFactor);				
		otherwise
			error('Wrongly Defined Line Direction!');
	end
	switch constantDir(2)
		case 'X'
			nodesOnLine = nodesOnLine1(find(abs(meshHierarchy_.boundaryNodeCoords(nodesOnLine1,1)-tarNode(1))<=lineRelaxFactor));
		case 'Y'
			nodesOnLine = nodesOnLine1(find(abs(meshHierarchy_.boundaryNodeCoords(nodesOnLine1,2)-tarNode(2))<=lineRelaxFactor));
		case 'Z'
			nodesOnLine = nodesOnLine1(find(abs(meshHierarchy_.boundaryNodeCoords(nodesOnLine1,3)-tarNode(3))<=lineRelaxFactor));				
		otherwise
			error('Wrongly Defined Line Direction!');
	end
	if 3==nargin
		endPoints = varargin{1}; endPoints = sort(endPoints, 'ascend');
		switch lineDir
			case 'X'
				nodesOnLine = nodesOnLine(meshHierarchy_.boundaryNodeCoords(nodesOnLine,1)>=endPoints(1));
				nodesOnLine = nodesOnLine(meshHierarchy_.boundaryNodeCoords(nodesOnLine,1)<=endPoints(2));
			case 'Y'
				nodesOnLine = nodesOnLine(meshHierarchy_.boundaryNodeCoords(nodesOnLine,2)>=endPoints(1));
				nodesOnLine = nodesOnLine(meshHierarchy_.boundaryNodeCoords(nodesOnLine,2)<=endPoints(2));					
			case 'Z'
				nodesOnLine = nodesOnLine(meshHierarchy_.boundaryNodeCoords(nodesOnLine,3)>=endPoints(1));
				nodesOnLine = nodesOnLine(meshHierarchy_.boundaryNodeCoords(nodesOnLine,3)<=endPoints(2));						
		end				
	end			

	hdpickedNode_(end+1) = plot3(axHandle, meshHierarchy_.boundaryNodeCoords(nodesOnLine,1), meshHierarchy_.boundaryNodeCoords(nodesOnLine,2), ...
		meshHierarchy_.boundaryNodeCoords(nodesOnLine,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);	
	numNewlyPickedNodes = length(nodesOnLine);
	pickedNodeCache_(end+1:end+numNewlyPickedNodes,1) = nodesOnLine;
end
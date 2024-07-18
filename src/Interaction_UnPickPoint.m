function Interaction_UnPickPoint(axHandle)
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdPickedNode_;
    dcm_obj = datacursormode;
    info_struct = getCursorInfo(dcm_obj);
    tarNode = info_struct.Position;
	if size(pickedNodeCache_,1)>0		
		potentialNodeList = meshHierarchy_(1).boundaryNodeCoords(pickedNodeCache_,:);
		[~, minValPos] = min(vecnorm(tarNode-potentialNodeList,2,2));	
		set(hdPickedNode_, 'visible', 'off');
		pickedNodeCache_(minValPos,:) = [];
		potentialNodeList(minValPos,:) = [];
		hdPickedNode_(end+1) = plot3(axHandle, potentialNodeList(:,1), potentialNodeList(:,2), ...
			potentialNodeList(:,3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);	
	end
end
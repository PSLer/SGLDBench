function Interaction_PickByPoint(axHandle, varargin)
	%% Syntax: 
	%% Interaction_PickByPoint(); 
	%% PickByLine(handInputCoordinate);
	global meshHierarchy_;
	global pickedNodeCache_;
	global hdPickedNode_;
	if 1==nargin
		dcm_obj = datacursormode;
		info_struct = getCursorInfo(dcm_obj);
		if isempty(info_struct)
			warning('No Cursor Mode Available!'); return;
		end
		tarNode = info_struct.Position;
	elseif 2==nargin
		tarNode = varargin{1}; tarNode = tarNode(:)';
	else
		error('Wrong Input!');		
	end
	hold(axHandle, 'on');
	hdPickedNode_(end+1) = plot3(axHandle, tarNode(1), tarNode(2), tarNode(3), 'xr', 'LineWidth', 2, 'MarkerSize', 6);
	[~,newlyPickedNode] = min(vecnorm(tarNode-meshHierarchy_(1).boundaryNodeCoords,2,2));	
	pickedNodeCache_(end+1,1) = newlyPickedNode;
end
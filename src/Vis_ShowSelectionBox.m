function Vis_ShowSelectionBox(cP1, cP2)
	global axHandle_;
	global hdSelectionBox_;
	x1 = cP1(1); y1 = cP1(2); z1 = cP1(3);
	x2 = cP2(1); y2 = cP2(2); z2 = cP2(3);
	
	selBox.vertices = [
		x1	y2	z1
		x1	y1	z1
		x2	y1	z1
		x2	y2	z1
		x1	y2	z2
		x1	y1	z2
		x2	y1	z2
		x2	y2	z2
	];
	selBox.faces = 1:8;
	selBox.faces = selBox.faces(:,[4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
	selBox.faces = reshape(selBox.faces, 4, 6)';
	set(hdSelectionBox_, 'visible', 'off');
	hold(axHandle_, 'on');
	hdSelectionBox_ = patch(axHandle_, selBox);
	set(hdSelectionBox_, 'faceColor', 'None', 'EdgeColor', 'magenta', 'lineWidth', 2);
	hold(axHandle_, 'on');
	hdSelectionBox_(end+1,1) = plot3(x1,y1,z1, '*y', 'MarkerSize', 10);
	hold(axHandle_, 'on');
	hdSelectionBox_(end+1,1) = plot3(x2,y2,z2, '+c', 'MarkerSize', 10);	
end
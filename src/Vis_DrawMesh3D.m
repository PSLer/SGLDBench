function hd = Vis_DrawMesh3D(axHandle, vertices, faces, edgeOpt)
	global meshHierarchy_;
	patchs2Draw.vertices = vertices;
	patchs2Draw.faces = faces;
	hd = patch(axHandle, patchs2Draw);
	if 0==edgeOpt
		set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1, 'EdgeColor', 'None');
		light(axHandle);
	else
		% if max([meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]) > 400
			% %% Omit mesh edges for high-resolution models
			% set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1, 'EdgeColor', 'None');
		% else
			% set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		% end		
		set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
	end	
	axis(axHandle, 'equal'); axis(axHandle, 'tight'); axis(axHandle, 'off');
end
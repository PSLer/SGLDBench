function Vis_ShowDesignByDensityLayoutInIsosurface(axHandle)
	global meshHierarchy_;
	global densityLayout_;

	passiveEles = find(densityLayout_>0.1);
	consideredEles = meshHierarchy_(1).eleMapBack(passiveEles);
	valForExtctIIsosurface = zeros(meshHierarchy_(1).resX*meshHierarchy_(1).resY*meshHierarchy_(1).resZ,1);			        
	valForExtctIIsosurface(consideredEles) = 1;
	valForExtctIIsosurface = reshape(valForExtctIIsosurface, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	valForExtctIIsosurface = flip(valForExtctIIsosurface,1);
	iSurface = isosurface(valForExtctIIsosurface,0);
	iCap = isocaps(valForExtctIIsosurface,0);
	hd = patch(axHandle, iSurface);
	hd(2) = patch(axHandle, iCap);
	set(hd, 'FaceColor', [255 240 214]/255, 'FaceAlpha', 1.0, 'EdgeColor', 'none');
	axis(axHandle, 'equal'); axis(axHandle, 'tight'); axis(axHandle, 'off');
end
function Vis_ShowDesignByDensityLayoutInIsosurface(axHandle)
	global meshHierarchy_;
	global densityLayout_;
	global voxelsOnBoundary_;
	global passiveElements_;
	
	disp('Extracting Isosurface of Design ...');
	solidEles = find(densityLayout_>0.1);
	solidElesExterior = [];
	if ~isempty(voxelsOnBoundary_)
		solidElesExterior = voxelsOnBoundary_;
	else
		if ~isempty(passiveElements_)		
			solidElesExterior = passiveElements_;
		end
	end
	solidEles = setdiff(solidEles, solidElesExterior);
	
	interiorEles = meshHierarchy_(1).eleMapBack(solidEles);
	valForExtctIIsosurface = zeros(meshHierarchy_(1).resX*meshHierarchy_(1).resY*meshHierarchy_(1).resZ,1);			        
	valForExtctIIsosurface(interiorEles) = 1;
	valForExtctIIsosurface = reshape(valForExtctIIsosurface, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	valForExtctIIsosurface = flip(valForExtctIIsosurface,1);
	iSurfaceInterior = isosurface(valForExtctIIsosurface,0);
	iCapInterior = isocaps(valForExtctIIsosurface,0);
	hd = patch(axHandle, iSurfaceInterior); hold(axHandle, 'on');
	hd(2) = patch(axHandle, iCapInterior);
	
	if ~isempty(solidElesExterior)	
		exteriorEles = meshHierarchy_(1).eleMapBack(solidElesExterior);
		valForExtctIIsosurface = zeros(meshHierarchy_(1).resX*meshHierarchy_(1).resY*meshHierarchy_(1).resZ,1);
		valForExtctIIsosurface(exteriorEles) = 1;
		valForExtctIIsosurface = reshape(valForExtctIIsosurface, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
		valForExtctIIsosurface = flip(valForExtctIIsosurface,1);
		iSurfaceExterior = isosurface(valForExtctIIsosurface,0);
		iCapExterior = isocaps(valForExtctIIsosurface,0);
		hold(axHandle, 'on'); hd1 = patch(axHandle, iSurfaceExterior);
		hold(axHandle, 'on'); hd1(2) = patch(axHandle, iCapExterior);
		set(hd1, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
	end
	set(hd, 'FaceColor', [255 240 214]/255, 'FaceAlpha', 1.0, 'EdgeColor', 'none');
	axis(axHandle, 'equal'); axis(axHandle, 'tight'); axis(axHandle, 'off');
	disp('Done with Isosurface Extraction!');
end
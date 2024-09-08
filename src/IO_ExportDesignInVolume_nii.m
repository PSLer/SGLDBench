function IO_ExportDesignInVolume_nii(fileName)
	global meshHierarchy_;
	global densityLayout_;
	global voxelsOnBoundary_;
	% V = zeros(numel(meshHierarchy_(1).eleMapForward),1, 'single');
	densityLayoutWithoutBoundary = densityLayout_;
	densityLayoutWithoutBoundary(voxelsOnBoundary_) = -1;
	V = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	V(meshHierarchy_(1).eleMapBack,1) = densityLayoutWithoutBoundary;
	V = reshape(V, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	% V = flip(V,1);
	niftiwrite(V,fileName);
end
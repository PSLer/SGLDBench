function IO_ExportDesignInVolume_nii(fileName)
	global meshHierarchy_;
	global densityLayoutWithoutBoundary_;
	V = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	V(meshHierarchy_(1).eleMapBack,1) = densityLayoutWithoutBoundary_;
	V = reshape(V, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	niftiwrite(V,fileName);
end
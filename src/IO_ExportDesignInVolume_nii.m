function IO_ExportDesignInVolume_nii(fileName)
	global meshHierarchy_;
	global densityLayout_;
	V = zeros(numel(meshHierarchy_(1).eleMapForward),1, 'single');
	V(meshHierarchy_(1).eleMapBack,1) = densityLayout_;
	V = reshape(V, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	niftiwrite(V,fileName);
end
function IO_ExportDesignInVolume_Geo_nii(fileName)
	global meshHierarchy_;
	global densityLayout4Vis_;
	V = reshape(densityLayout4Vis_, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	% V = flip(V,1);
	niftiwrite(V,fileName);
end
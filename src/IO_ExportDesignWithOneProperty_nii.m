function IO_ExportDesignWithOneProperty_nii(propertyVolume, fileName)
	global meshHierarchy_;
	global densityLayout_;
	global densityLayout4Vis_;
	
	if isempty(densityLayout4Vis_)
		Vdensity = zeros(numel(meshHierarchy_(1).eleMapForward),1);
		Vdensity(meshHierarchy_(1).eleMapBack,1) = densityLayout_;
		Vdensity = reshape(Vdensity, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);		
	else
		Vdensity = densityLayout4Vis_;
		Vdensity = reshape(Vdensity, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);	
	end
	V = zeros(meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ, 2);
	V(:,:,:,1) = Vdensity;
	V(:,:,:,2) = propertyVolume;
	V = permute(V, [4, 1, 2, 3]);
	
	% V = flip(V,1);
	niftiwrite(V,fileName);
end

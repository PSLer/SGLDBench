function IO_ExportDesignWithOneProperty_nii(propertyVolume, fileName, varargin)
	global meshHierarchy_;
	global densityLayout_;
	global voxelsOnBoundary_;
	densityLayoutWithoutBoundary = densityLayout_;
	if 3==nargin
		if varargin{1}
			densityLayoutWithoutBoundary(voxelsOnBoundary_) = -1;
		end
	end
	Vdensity = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	Vdensity(meshHierarchy_(1).eleMapBack,1) = densityLayoutWithoutBoundary;
	Vdensity = reshape(Vdensity, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	
	V = zeros(meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ, 2);
	V(:,:,:,1) = Vdensity;
	V(:,:,:,2) = propertyVolume;
	V = permute(V, [4, 1, 2, 3]);
	
	% V = flip(V,1);
	niftiwrite(V,fileName);
end

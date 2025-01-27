function perEleVolume = Common_ConvertPerEleVector2Volume(perEleVec)
	global meshHierarchy_;
	perEleVolume = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	perEleVolume(meshHierarchy_(1).eleMapBack,1) = perEleVec;
	perEleVolume = reshape(perEleVolume, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
end
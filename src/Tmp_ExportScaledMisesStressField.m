function Tmp_ExportScaledMisesStressField()
	global meshHierarchy_;
	global cartesianStressField_;
	
	stressAtCellCtrs = zeros(meshHierarchy_(1).numElements,6);
	shapeFuncsAtCentroid = ones(1,8)/8;
	for ii=1:meshHierarchy_(1).numElements
		iStress = cartesianStressField_(meshHierarchy_(1).eNodMat(ii,:),:);
		stressAtCellCtrs(ii,:) = shapeFuncsAtCentroid * iStress;
	end
	
	misesStressAtCellCtrs = sqrt(0.5*((stressAtCellCtrs(:,1)-stressAtCellCtrs(:,2)).^2 + ...
		(stressAtCellCtrs(:,2)-stressAtCellCtrs(:,3)).^2 + (stressAtCellCtrs(:,3)...
			-stressAtCellCtrs(:,1)).^2 ) + 3*( stressAtCellCtrs(:,6).^2 + stressAtCellCtrs(:,4).^2 + ...
				stressAtCellCtrs(:,5).^2 ));
				
misesStressAtCellCtrs = misesStressAtCellCtrs.^(1/1);

	misesStressVolume = zeros(numel(meshHierarchy_(1).eleMapForward),1);
	misesStressVolume(meshHierarchy_(1).eleMapBack,1) = misesStressAtCellCtrs;
	misesStressVolume = reshape(misesStressVolume, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	niftiwrite(misesStressVolume, './out/vonMisesStressVolume.nii');
end
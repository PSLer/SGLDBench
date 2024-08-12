function c = CmptComplianceTemp()
	global meshHierarchy_;
    global U_;
	blockIndex = Solving_MissionPartition(meshHierarchy_(1).numElements, 1.0e7);
	ceList = zeros(meshHierarchy_(1).numElements, 1);
	for ii=1:size(blockIndex,1)
		rangeIndex = (blockIndex(ii,1):blockIndex(ii,2))';
		iReshapedU = zeros(numel(rangeIndex),24);
		% iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
		% iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
        iElesNodMat = meshHierarchy_(1).eNodMat(rangeIndex,:);
        tmp = U_(1:3:end,:); iReshapedU(:,1:3:24) = tmp(iElesNodMat);
        tmp = U_(2:3:end,:); iReshapedU(:,2:3:24) = tmp(iElesNodMat);
        tmp = U_(3:3:end,:); iReshapedU(:,3:3:24) = tmp(iElesNodMat);		
		ceList(rangeIndex,1) = sum((iReshapedU*meshHierarchy_(1).Ke).*iReshapedU,2);
	end	
    c = meshHierarchy_(1).eleModulus(:) .* ceList(:); c =sum(c);
end
function ceList = TopOpti_ComputeUnitCompliance(printLSS_Proc)
	global meshHierarchy_;
	global U_;
	
	Solving_AssembleFEAstencil();
	Solving_CG_GMGS(printLSS_Proc);
	
	blockIndex = MissionPartition(meshHierarchy_(1).numElements, 1.0e6);
	ceList = zeros(meshHierarchy_(1).numElements, 1);
	for ii=1:size(blockIndex,1)
		rangeIndex = (blockIndex(ii,1):blockIndex(ii,2))';
		iReshapedU = zeros(numel(rangeIndex),24);
        tmp = U_(1:3:end,:); iReshapedU(:,1:3:24) = tmp(meshHierarchy_(1).eNodMat(rangeIndex,:));
        tmp = U_(2:3:end,:); iReshapedU(:,2:3:24) = tmp(meshHierarchy_(1).eNodMat(rangeIndex,:));
        tmp = U_(3:3:end,:); iReshapedU(:,3:3:24) = tmp(meshHierarchy_(1).eNodMat(rangeIndex,:));		
		ceList(rangeIndex,1) = sum((iReshapedU*meshHierarchy_(1).Ke).*iReshapedU,2);
	end		
end
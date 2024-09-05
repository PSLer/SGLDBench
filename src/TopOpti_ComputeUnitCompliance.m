function ceList = TopOpti_ComputeUnitCompliance()
	global MEXfunc_;
    global meshHierarchy_;
	global U_;	
	
    if MEXfunc_
        ceList = TopOpti_CmptUnitCompliance_mex(U_, meshHierarchy_(1).eNodMat, meshHierarchy_(1).Ke, 1.0e7);
    else
	    blockIndex = Solving_MissionPartition(meshHierarchy_(1).numElements, 1.0e7);
	    ceList = zeros(meshHierarchy_(1).numElements, 1);
	    for ii=1:size(blockIndex,1)
		    rangeIndex = (blockIndex(ii,1):blockIndex(ii,2))';
		    iReshapedU = zeros(numel(rangeIndex),24);
		    iElesNodMat = meshHierarchy_(1).eNodMat(rangeIndex,:);
            tmp = U_(1:3:end,:); iReshapedU(:,1:3:24) = tmp(iElesNodMat);
            tmp = U_(2:3:end,:); iReshapedU(:,2:3:24) = tmp(iElesNodMat);
            tmp = U_(3:3:end,:); iReshapedU(:,3:3:24) = tmp(iElesNodMat);		
		    ceList(rangeIndex,1) = sum((iReshapedU*meshHierarchy_(1).Ke).*iReshapedU,2);
	    end	
    end	
end
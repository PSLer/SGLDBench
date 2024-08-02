function rCoaser = Solving_RestrictResidual(rFiner,ii)
	global meshHierarchy_;
	rFiner = reshape(rFiner,3,meshHierarchy_(ii-1).numNodes)';
	rCoaser = zeros(meshHierarchy_(ii).numNodes,3);
	rFiner1 = zeros(meshHierarchy_(ii).intermediateNumNodes,3);
	rFiner1(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:) = rFiner;
	rFiner1 = rFiner1./meshHierarchy_(ii).transferMatCoeffi;
	eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(ii).eNodMatHalf);
	eNodMat = eNodMat(:);
	for jj=1:3
		tmp = rFiner1(:,jj);
		tmp = tmp(meshHierarchy_(ii).transferMat);
		tmp = tmp' * meshHierarchy_(ii).multiGridOperatorRI;
		rCoaser(:,jj) = accumarray(eNodMat,tmp(:),[meshHierarchy_(ii).numNodes 1]);
	end	
	rCoaser = reshape(rCoaser', 3*meshHierarchy_(ii).numNodes, 1);	
end
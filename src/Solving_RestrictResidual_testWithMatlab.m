function rCoaser = Solving_RestrictResidual_testWithMatlab(rFiner,ii)
	global meshHierarchy_;
	global MEXfunc_;
dimN = size(rFiner,1) / meshHierarchy_(ii-1).numNodes;	
	rFiner = reshape(rFiner,dimN,meshHierarchy_(ii-1).numNodes)';
	rFiner1 = zeros(meshHierarchy_(ii).intermediateNumNodes,dimN);
	rFiner1(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:) = rFiner;
	rFiner1 = rFiner1./meshHierarchy_(ii).transferMatCoeffi;
	for jj=1:dimN
		tmp = rFiner1(:,jj);
		tmp = tmp(meshHierarchy_(ii).transferMat);
		tmp = tmp' * meshHierarchy_(ii).multiGridOperatorRI;
		rCoaser(:,jj) = accumarray(meshHierarchy_(ii).eNodMat(:),tmp(:),[meshHierarchy_(ii).numNodes 1]);
	end	
	rCoaser = reshape(rCoaser', dimN*meshHierarchy_(ii).numNodes, 1);
end
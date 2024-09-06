function rCoaser = Solving_RestrictResidual(rFiner,ii)
	global meshHierarchy_;
	global MEXfunc_;
	rFiner = reshape(rFiner,3,meshHierarchy_(ii-1).numNodes)';
	if MEXfunc_
		rCoaser = Solving_Restriction_MatrixFree_mex_advanced(rFiner, meshHierarchy_(ii).transferMat, meshHierarchy_(ii).multiGridOperatorRI, ...
			meshHierarchy_(ii).eNodMat, meshHierarchy_(ii).numNodes, meshHierarchy_(ii).intermediateNumNodes, meshHierarchy_(ii).solidNodeMapCoarser2Finer, meshHierarchy_(ii).transferMatCoeffi);	
	else
		rFiner1 = zeros(meshHierarchy_(ii).intermediateNumNodes,3);
		rFiner1(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:) = rFiner;
		rFiner1 = rFiner1./meshHierarchy_(ii).transferMatCoeffi;
		rCoaser = zeros(meshHierarchy_(ii).numNodes,3);		
		for jj=1:3
			tmp = rFiner1(:,jj);
			tmp = tmp(meshHierarchy_(ii).transferMat);
			tmp = tmp' * meshHierarchy_(ii).multiGridOperatorRI;
			rCoaser(:,jj) = accumarray(meshHierarchy_(ii).eNodMat(:),tmp(:),[meshHierarchy_(ii).numNodes 1]);
		end	
		clear rFiner1
	end
	rCoaser = reshape(rCoaser', 3*meshHierarchy_(ii).numNodes, 1);	
end
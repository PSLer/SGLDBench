function rCoaser = Solving_RestrictResidual(rFiner,ii)
	global meshHierarchy_;
	global MEXfunc_;
	rFiner = reshape(rFiner,3,meshHierarchy_(ii-1).numNodes)';
	rFiner1 = zeros(meshHierarchy_(ii).intermediateNumNodes,3);
	rFiner1(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:) = rFiner;
	rFiner1 = rFiner1./meshHierarchy_(ii).transferMatCoeffi;
	if 0
		rCoaser = Solving_Restriction_MatrixFree_mex(rFiner1, meshHierarchy_(ii).transferMat', meshHierarchy_(ii).multiGridOperatorRI, meshHierarchy_(ii).eNodMat, meshHierarchy_(ii).numNodes);
	else
		rCoaser = zeros(meshHierarchy_(ii).numNodes,3);
		% eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(ii).eNodMatHalf);
		% eNodMat = eNodMat(:);
		mpOpt = 1; %%1==Previous, 0==New
		if mpOpt
			for jj=1:3
				tmp = rFiner1(:,jj);
				tmp = tmp(meshHierarchy_(ii).transferMat);
				tmp = tmp' * meshHierarchy_(ii).multiGridOperatorRI;
				rCoaser(:,jj) = accumarray(meshHierarchy_(ii).eNodMat(:),tmp(:),[meshHierarchy_(ii).numNodes 1]);
			end	
		else
			for jj=1:3
				tmp = rFiner1(:,jj);
				tmp = Vector2Matrix_Indexing_mex(tmp, meshHierarchy_(ii).transferMat');
				tmp = tmp * meshHierarchy_(ii).multiGridOperatorRI;
				rCoaser(:,jj) = Accumarray_mex(meshHierarchy_(ii).eNodMat(:),tmp(:),[meshHierarchy_(ii).numNodes 1]);
			end	
		end		
	end
	rCoaser = reshape(rCoaser', 3*meshHierarchy_(ii).numNodes, 1);
end
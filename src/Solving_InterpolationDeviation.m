function xFiner = Solving_InterpolationDeviation(xCoarser, ii)
	global meshHierarchy_;
	global MEXfunc_;
	xCoarser = reshape(xCoarser,3,meshHierarchy_(ii).numNodes)';
	xFiner = zeros(meshHierarchy_(ii).intermediateNumNodes,3);
	eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(ii).eNodMatHalf);
	transferMat = meshHierarchy_(ii).transferMat(:);
	impOpt = 0; %%1==Previous, 0==New
	if MEXfunc_
		for jj=1:3
			tmp = xCoarser(:,jj);
			tmp = Vector2Matrix_Indexing_mex(tmp, eNodMat); 
			tmp1 = meshHierarchy_(ii).multiGridOperatorRI * tmp';
			xFiner(:,jj) = Accumarray_mex(transferMat,tmp1(:),[meshHierarchy_(ii).intermediateNumNodes 1]);
		end			
	else
		for jj=1:3
			tmp = xCoarser(:,jj);
			tmp = tmp(eNodMat);
			tmp1 = meshHierarchy_(ii).multiGridOperatorRI * tmp';
			xFiner(:,jj) = accumarray(transferMat,tmp1(:),[meshHierarchy_(ii).intermediateNumNodes 1]);
		end
	end
	xFiner = xFiner ./ meshHierarchy_(ii).transferMatCoeffi;
	xFiner = xFiner(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:);
	xFiner = reshape(xFiner', 3*meshHierarchy_(ii-1).numNodes, 1);
end
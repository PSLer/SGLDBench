function xFiner = Solving_InterpolationDeviation(xCoarser, ii)
	global meshHierarchy_;
	global MEXfunc_;
	
	xCoarser = reshape(xCoarser,3,meshHierarchy_(ii).numNodes)';
	if MEXfunc_
		xFiner = Solving_Interpolation_MatrixFree_mex(xCoarser, meshHierarchy_(ii).transferMat, meshHierarchy_(ii).multiGridOperatorRI, ...
			meshHierarchy_(ii).eNodMat, meshHierarchy_(ii-1).numNodes, meshHierarchy_(ii).intermediateNumNodes, meshHierarchy_(ii).solidNodeMapCoarser2Finer, meshHierarchy_(ii).transferMatCoeffi);		
	else
		
		xFiner = zeros(meshHierarchy_(ii).intermediateNumNodes,3);
		transferMat = meshHierarchy_(ii).transferMat(:);
		impOpt = 0; %%True==IndexingMEX, false==pureMatlab
		if impOpt
			for jj=1:3
				tmp = xCoarser(:,jj);
				tmp = Vector2Matrix_Indexing_mex(tmp, meshHierarchy_(ii).eNodMat); 
				tmp1 = meshHierarchy_(ii).multiGridOperatorRI * tmp';
				xFiner(:,jj) = Accumarray_mex(transferMat,tmp1(:),[meshHierarchy_(ii).intermediateNumNodes 1]);
			end			
		else
			for jj=1:3
				tmp = xCoarser(:,jj);
				tmp = tmp(meshHierarchy_(ii).eNodMat);
				tmp1 = meshHierarchy_(ii).multiGridOperatorRI * tmp';
				xFiner(:,jj) = accumarray(transferMat,tmp1(:),[meshHierarchy_(ii).intermediateNumNodes 1]);
			end
		end	
		xFiner = xFiner ./ meshHierarchy_(ii).transferMatCoeffi;
		xFiner = xFiner(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:);		
	end
	xFiner = reshape(xFiner', 3*meshHierarchy_(ii-1).numNodes, 1);
end
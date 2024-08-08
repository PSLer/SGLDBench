function rCoaser = Solving_RestrictResidual_testWithMEX(rFiner,ii)
	global meshHierarchy_;
	global MEXfunc_;
	rFiner = reshape(rFiner,3,meshHierarchy_(ii-1).numNodes)';
	rFiner1 = zeros(meshHierarchy_(ii).intermediateNumNodes,3);
	rFiner1(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:) = rFiner;
	rFiner1 = rFiner1./meshHierarchy_(ii).transferMatCoeffi;
	transferMatT = meshHierarchy_(ii).transferMat';
	for jj=1:3
		rCoaser(:,jj) = Solving_Restriction_MatrixFree_mex(rFiner1(:,jj), transferMatT, meshHierarchy_(ii).multiGridOperatorRI, meshHierarchy_(ii).eNodMat, meshHierarchy_(ii).numNodes);
	end
	rCoaser = reshape(rCoaser', 3*meshHierarchy_(ii).numNodes, 1);
end
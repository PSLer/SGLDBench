function rCoaser = Solving_RestrictResidual_testWithMatlab(rFiner1,ii)
	global meshHierarchy_;
if 1 %%test
	rCoaser = rFiner1(meshHierarchy_(ii).transferMat);
	rCoaser = rCoaser' * meshHierarchy_(ii).multiGridOperatorRI;
	rCoaser = accumarray(meshHierarchy_(ii).eNodMat(:),rCoaser(:),[meshHierarchy_(ii).numNodes 1]);
else	
	tmp = rFiner1(meshHierarchy_(ii).transferMat);
	tmp = tmp' * meshHierarchy_(ii).multiGridOperatorRI;
	rCoaser = accumarray(meshHierarchy_(ii).eNodMat(:),tmp(:),[meshHierarchy_(ii).numNodes 1]);
end
end
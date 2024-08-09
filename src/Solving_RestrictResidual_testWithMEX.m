function rCoaser = Solving_RestrictResidual_testWithMEX(rFiner,ii)
	global meshHierarchy_;
	global MEXfunc_;
	rFiner = reshape(rFiner,3,meshHierarchy_(ii-1).numNodes)';
	if 1
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
	end
	rCoaser = reshape(rCoaser', 3*meshHierarchy_(ii).numNodes, 1);	
end

function rCoaser = Solving_RestrictResidual(rFiner,ii)
	global meshHierarchy_;
	global MEXfunc_;
tStart1 = tic;	

tEnd1 = toc(tStart1);	
tStart2 = tic;	
	if 1
		rCoaser = Solving_Restriction_MatrixFree_mex(rFiner1, meshHierarchy_(ii).transferMat, meshHierarchy_(ii).multiGridOperatorRI, meshHierarchy_(ii).eNodMat, meshHierarchy_(ii).numNodes);			
	else
		mpOpt = 1; %%1==Previous, 0==New
		if mpOpt

		else
			for jj=1:3
				tmp = rFiner1(:,jj);
				tmp = Vector2Matrix_Indexing_mex(tmp, meshHierarchy_(ii).transferMat');
				tmp = tmp * meshHierarchy_(ii).multiGridOperatorRI;
				rCoaser(:,jj) = Accumarray_mex(meshHierarchy_(ii).eNodMat(:),tmp(:),[meshHierarchy_(ii).numNodes 1]);
			end	
		end
	end
tEnd2 = toc(tStart2);
tStart3 = tic;

tEnd3 = toc(tStart3);
disp(['Data Preparation Costs:' sprintf('%f', tEnd1)]);
disp(['Perform Computation Costs:' sprintf('%f', tEnd2)]);
disp(['Finalize Output Costs:' sprintf('%f', tEnd3)]);
disp(['Total Costs:' sprintf('%f', toc(tStart1))]);
end
function xFiner = Solving_InterpolationDeviation(xCoarser, ii)
	global meshHierarchy_;
	
	xCoarser = reshape(xCoarser,3,meshHierarchy_(ii).numNodes)';
	xFiner = zeros(meshHierarchy_(ii).intermediateNumNodes,3);
	nn = (meshHierarchy_(ii).spanWidth+1)^3;
	impOpt = 0; %%1==Previous, 0==New
	if impOpt
		for jj=1:3
			tmp = xCoarser(:,jj);
			tmp = tmp(Common_RecoverHalfeNodMat(meshHierarchy_(ii).eNodMatHalf));
			tmp1 = meshHierarchy_(ii).multiGridOperatorRI * tmp';
			for kk=1:nn
				xFiner(meshHierarchy_(ii).transferMat(kk,:),jj) = ...
					xFiner(meshHierarchy_(ii).transferMat(kk,:),jj) + tmp1(kk,:)';
			end			
		end	
	else
		for jj=1:3
			tmp = xCoarser(:,jj);
			tmp = tmp(Common_RecoverHalfeNodMat(meshHierarchy_(ii).eNodMatHalf));
			tmp1 = meshHierarchy_(ii).multiGridOperatorRI * tmp';
			xFiner(:,jj) = accumarray(meshHierarchy_(ii).transferMat(:),tmp1(:),[meshHierarchy_(ii).intermediateNumNodes 1]);
		end	
	end

	xFiner = xFiner ./ meshHierarchy_(ii).transferMatCoeffi;
	xFiner = xFiner(meshHierarchy_(ii).solidNodeMapCoarser2Finer,:);
	xFiner = reshape(xFiner', 3*meshHierarchy_(ii-1).numNodes, 1);
end
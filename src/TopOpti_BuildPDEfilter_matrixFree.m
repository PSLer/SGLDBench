function TopOpti_BuildPDEfilter_matrixFree()
	global meshHierarchy_;
	global rHatMin_; 
	global KePDE_;
	global diagKePDE_;
	
	[s, t, p, w] = FEA_GaussianIntegral();
	N = FEA_ShapeFunction(s, t, p);
	dShape = FEA_DeShapeFunction(s,t,p);
	
	%%Jacobian Matrix, corresponding to the commonly used 2x2x2 cubic element in natural coordinate system
	%%CellSize
	CellSize = 1; %%Always treated as a unit cell
	detJ = CellSize^3 /8 * ones(8,1); %%Sub-Volume
	wgt = w.*detJ;
	KEF0 = dShape'*dShape;
	KEF1 = N'*diag(wgt)*N;
    iRmin = (rHatMin_ * meshHierarchy_(1).eleSize(1))/2/sqrt(3);
    iKEF = iRmin^2*KEF0 + KEF1;
	KePDE_ = iKEF;
	
	%%Diagonal Preconditioner
	diagKePDE_ = zeros(meshHierarchy_(1).numNodes,1);
	numElements = meshHierarchy_(1).numElements;
	diagKe = diag(KePDE_);
	blockIndex = Solving_MissionPartition(numElements, 1.0e7);
	for jj=1:size(blockIndex,1)
		rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))';
		% jElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
		% jElesNodMat = Common_RecoverHalfeNodMat(jElesNodMat)';
		jElesNodMat = meshHierarchy_(1).eNodMat(rangeIndex,:)';
		diagKeBlock = diagKe(:) .* ones(1,numel(rangeIndex));
		jElesNodMat = jElesNodMat(:);
		diagKeBlockSingleDOF = diagKeBlock(:); 
		diagKePDE_ = diagKePDE_ + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(1).numNodes, 1]);				
	end	
	diagKePDE_ = diagKePDE_.^(-1);
	
	% %%Assembling Exclusive Computing Stencil
	% supMeshHierarchy4PDEfilter_(1).Ke = iKEF;
	% supMeshHierarchy4PDEfilter_(1).Ks = supMeshHierarchy4PDEfilter_(1).Ke;
	% for ii=2:numel(meshHierarchy_)
		% spanWidth = meshHierarchy_(ii).spanWidth;
		% interpolatingKe = Solving_Operator4MultiGridRestrictionAndInterpolation('inNODE',spanWidth);
		% eNodMat4Finer2Coarser = Solving_SubEleNodMat(spanWidth);
		% [rowIndice, colIndice, ~] = find(ones(8));
		% iK = eNodMat4Finer2Coarser(:,rowIndice)';
		% jK = eNodMat4Finer2Coarser(:,colIndice)';
		% numProjectDOFs = (spanWidth+1)^3; %% 1 DOF == 1 NODE
		% supMeshHierarchy4PDEfilter_(ii).storingState = 1;
		% supMeshHierarchy4PDEfilter_(ii).Ke = supMeshHierarchy4PDEfilter_(ii-1).Ke*spanWidth;
		% numElements = meshHierarchy_(ii).numElements;
		% Ks = repmat(supMeshHierarchy4PDEfilter_(ii).Ke, 1,1,numElements);
		% finerKes = zeros(8*8,spanWidth^3);
		% elementUpwardMap = meshHierarchy_(ii).elementUpwardMap;
		% if 2==ii
			% iKe = supMeshHierarchy4PDEfilter_(ii-1).Ke;
			% iKs = reshape(iKe, 8*8, 1);
			% parpool('threads');	
			% parfor jj=1:numElements
				% sonEles = elementUpwardMap(jj,:);
				% solidEles = find(0~=sonEles);
				% sK = finerKes;
				% sK(:,solidEles) = repmat(iKs, 1, numel(solidEles));
				% tmpK = sparse(iK, jK, sK, numProjectDOFs, numProjectDOFs);
				% tmpK = interpolatingKe' * tmpK * interpolatingKe;
				% Ks(:,:,jj) = full(tmpK);
			% end
			% delete(gcp('nocreate')); 
		% else
			% KsPrevious = supMeshHierarchy4PDEfilter_(ii-1).Ks;
			% parpool('threads');
			% parfor jj=1:numElements
				% iFinerEles = elementUpwardMap(jj,:);
				% solidEles = find(0~=iFinerEles);
				% iFinerEles = iFinerEles(solidEles);
				% sK = finerKes;
				% tarKes = KsPrevious(:,:,iFinerEles);
				% for kk=1:length(solidEles)
					% sK(:,solidEles(kk)) = reshape(tarKes(:,:,kk),8^2,1);
				% end
				% tmpK = sparse(iK, jK, sK, numProjectDOFs, numProjectDOFs);
				% tmpK = interpolatingKe' * tmpK * interpolatingKe;
				% Ks(:,:,jj) = full(tmpK);
			% end			
			% delete(gcp('nocreate'));
		% end
		% supMeshHierarchy4PDEfilter_(ii).Ks = Ks;
	% end
	
	% %% initialize smoother
	% for ii=1:length(meshHierarchy_)-1
		% diagK = zeros(meshHierarchy_(ii).numNodes,1);
		% numElements = meshHierarchy_(ii).numElements;
		% Ks = supMeshHierarchy4PDEfilter_(ii).Ks;
		% if 1==size(Ks,3)
			% diagKe = diag(supMeshHierarchy4PDEfilter_(ii).Ks);
			% blockIndex = Solving_MissionPartition(numElements, 1.0e7);
			% for jj=1:size(blockIndex,1)
				% rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))';
				% jElesNodMat = meshHierarchy_(ii).eNodMatHalf(rangeIndex,:);
				% jElesNodMat = Common_RecoverHalfeNodMat(jElesNodMat)';
				% diagKeBlock = diagKe(:) .* ones(1,numel(rangeIndex));
				% jElesNodMat = jElesNodMat(:);
				% diagKeBlockSingleDOF = diagKeBlock(:); 
				% diagK = diagK + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(ii).numNodes, 1]);				
			% end
		% else
			% blockIndex = Solving_MissionPartition(numElements, 1.0e7);
			% for jj=1:size(blockIndex,1)
				% rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))';
				% jElesNodMat = meshHierarchy_(ii).eNodMatHalf(rangeIndex,:);
				% jElesNodMat = Common_RecoverHalfeNodMat(jElesNodMat)';
				% jKs = Ks(:,:,rangeIndex);
				% jKs = reshape(jKs,8*8,numel(rangeIndex));
				% diagKeBlock = jKs(1:9:(8*8),:);
				% jElesNodMat = jElesNodMat(:);
				% diagKeBlockSingleDOF = diagKeBlock(:);
				% diagK = diagK + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(ii).numNodes, 1]);
			% end
		% end
		% supMeshHierarchy4PDEfilter_(ii).diagK = diagK;
	% end
	
	% %% Assemble&Factorize PDF Kernel Matrix on Coarsest Level
	% [rowIndice, colIndice, ~] = find(ones(8));	
	% sK = zeros(8^2, meshHierarchy_(end).numElements);
	% for ii=1:meshHierarchy_(end).numElements
		% sK(:,ii) = reshape(supMeshHierarchy4PDEfilter_(end).Ks(:,:,ii), 8^2, 1);
	% end
	% eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(end).eNodMatHalf);
	% iK = eNodMat(:,rowIndice);
	% jK = eNodMat(:,colIndice);
	% KcoarsestLevel = sparse(iK, jK, sK');
	% [cholFacPDE_, ~, cholPermutPDE_] = chol(KcoarsestLevel,'lower');
end
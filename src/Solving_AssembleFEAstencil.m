function Solving_AssembleFEAstencil()
	global meshHierarchy_;
	global numLevels_;
	global cholFac_; global cholPermut_;
	global MEXfunc_;
	
	%% Compute 'Ks' on Coarser Levels
	reOrdering = [1 9 17 2 10 18 3 11 19 4 12 20 5 13 21 6 14 22 7 15 23 8 16 24];
	for ii=2:numLevels_
		spanWidth = meshHierarchy_(ii).spanWidth;
		interpolatingKe = Solving_Operator4MultiGridRestrictionAndInterpolation('inDOF',spanWidth);
		eNodMat4Finer2Coarser = Solving_SubEleNodMat(spanWidth);
		[rowIndice, colIndice, ~] = find(ones(24));
		eDofMat4Finer2Coarser = [3*eNodMat4Finer2Coarser-2 3*eNodMat4Finer2Coarser-1 3*eNodMat4Finer2Coarser];
		eDofMat4Finer2Coarser = eDofMat4Finer2Coarser(:, reOrdering);
		iK = eDofMat4Finer2Coarser(:,rowIndice)';
		jK = eDofMat4Finer2Coarser(:,colIndice)';
		numProjectNodes = (spanWidth+1)^3;
		numProjectDOFs = numProjectNodes*3;
		localMapping = iK(:) + (jK(:)-1)*numProjectDOFs; localMapping = int32(localMapping);
		meshHierarchy_(ii).storingState = 1;
		meshHierarchy_(ii).Ke = meshHierarchy_(ii-1).Ke*spanWidth;
		numElements = meshHierarchy_(ii).numElements;	
		diagK = zeros(meshHierarchy_(ii).numNodes,3);
		finerKes = zeros(24*24,spanWidth^3);
		elementUpwardMap = meshHierarchy_(ii).elementUpwardMap;
		%%Compute Element Stiffness Matrices on Coarser Levels
		if 2==ii			
			iKe = meshHierarchy_(ii-1).Ke;
			iKs = reshape(iKe, 24*24, 1);
			eleModulus = meshHierarchy_(1).eleModulus;
			if MEXfunc_							
				Ks = Solving_AssembleCmptStencilFromFinestLevel(iKe, eleModulus, elementUpwardMap, interpolatingKe, localMapping, numProjectNodes);			
			else
				Ks = repmat(meshHierarchy_(ii).Ke, 1,1,numElements);
				if isempty(gcp('nocreate')), parpool('threads'); end			
				parfor jj=1:numElements
					sonEles = elementUpwardMap(jj,:);
					solidEles = find(0~=sonEles);
					sK = finerKes;
					sK(:,solidEles) = iKs .* eleModulus(sonEles(solidEles));
					%%previous				
					tmpK = sparse(iK, jK, sK, numProjectDOFs, numProjectDOFs);
					tmpK = interpolatingKe' * tmpK * interpolatingKe;
					Ks(:,:,jj) = full(tmpK);				
					%%New slightly faster
					% tmpK = accumarray(localMapping, sK(:), [numProjectDOFs^2, 1]); 
					% tmpK = reshape(tmpK, numProjectDOFs, numProjectDOFs);
					% Ks(:,:,jj) = interpolatingKe' * tmpK * interpolatingKe;
				end			
			end
		else
			% KsPrevious = meshHierarchy_(ii-1).Ks;
            KsPrevious = Ks; clear Ks;
			if MEXfunc_			
				Ks = Solving_AssembleCmptStencilFromNonFinestLevel(KsPrevious, elementUpwardMap, interpolatingKe, localMapping, numProjectNodes);								
			else
				Ks = repmat(meshHierarchy_(ii).Ke, 1,1,numElements);
				if isempty(gcp('nocreate')), parpool('threads'); end	
				parfor jj=1:numElements
					iFinerEles = elementUpwardMap(jj,:);
					solidEles = find(0~=iFinerEles);
					iFinerEles = iFinerEles(solidEles);
					sK = finerKes;
					tarKes = KsPrevious(:,:,iFinerEles);
					for kk=1:length(solidEles)
						sK(:,solidEles(kk)) = reshape(tarKes(:,:,kk),24^2,1);
					end
					%%previous
					tmpK = sparse(iK, jK, sK, numProjectDOFs, numProjectDOFs);
					tmpK = interpolatingKe' * tmpK * interpolatingKe;
					Ks(:,:,jj) = full(tmpK);				
					%%New
					% tmpK = accumarray(localMapping, sK(:), [numProjectDOFs^2, 1]); 
					% tmpK = reshape(tmpK, numProjectDOFs, numProjectDOFs);
					% Ks(:,:,jj) = interpolatingKe' * tmpK * interpolatingKe;
				end			
			end				
		end			
		% meshHierarchy_(ii).Ks = Ks;
		%%Initialize Jacobian Smoother on Coarser Levels
		if ii<numLevels_
			blockIndex = Solving_MissionPartition(numElements, 1.0e7);
			for jj=1:size(blockIndex,1)
				rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))';
				jElesNodMat = meshHierarchy_(ii).eNodMat(rangeIndex,:)';
				jKs = Ks(:,:,rangeIndex);
				jKs = reshape(jKs,24*24,numel(rangeIndex));
				diagKeBlock = jKs(1:25:(24*24),:);
				jElesNodMat = jElesNodMat(:);
				diagKeBlockSingleDOF = diagKeBlock(1:3:end,:); diagKeBlockSingleDOF = diagKeBlockSingleDOF(:);
				diagK(:,1) = diagK(:,1) + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(ii).numNodes, 1]);
				diagKeBlockSingleDOF = diagKeBlock(2:3:end,:); diagKeBlockSingleDOF = diagKeBlockSingleDOF(:);
				diagK(:,2) = diagK(:,2) + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(ii).numNodes, 1]);
				diagKeBlockSingleDOF = diagKeBlock(3:3:end,:); diagKeBlockSingleDOF = diagKeBlockSingleDOF(:);
				diagK(:,3) = diagK(:,3) + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(ii).numNodes, 1]);		
			end
			meshHierarchy_(ii).diagK = reshape(diagK',meshHierarchy_(ii).numDOFs,1);
        end		
	end		

	%%%Initialize Jacobian Smoother on Finest Level
	diagK = zeros(meshHierarchy_(1).numNodes,3);
	numElements = meshHierarchy_(1).numElements;
	diagKe = diag(meshHierarchy_(1).Ke);
	eleModulus = meshHierarchy_(1).eleModulus;
	blockIndex = Solving_MissionPartition(numElements, 1.0e7);
	for jj=1:size(blockIndex,1)				
		rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))';
		jElesNodMat = meshHierarchy_(1).eNodMat(rangeIndex,:)';
		jEleModulus = eleModulus(1, rangeIndex);
		diagKeBlock = diagKe(:) .* jEleModulus;
		jElesNodMat = jElesNodMat(:);
		diagKeBlockSingleDOF = diagKeBlock(1:3:end,:); diagKeBlockSingleDOF = diagKeBlockSingleDOF(:);
		diagK(:,1) = diagK(:,1) + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(1).numNodes, 1]);
		diagKeBlockSingleDOF = diagKeBlock(2:3:end,:); diagKeBlockSingleDOF = diagKeBlockSingleDOF(:);
		diagK(:,2) = diagK(:,2) + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(1).numNodes, 1]);
		diagKeBlockSingleDOF = diagKeBlock(3:3:end,:); diagKeBlockSingleDOF = diagKeBlockSingleDOF(:);
		diagK(:,3) = diagK(:,3) + accumarray(jElesNodMat, diagKeBlockSingleDOF, [meshHierarchy_(1).numNodes, 1]);	
	end
	meshHierarchy_(1).diagK = reshape(diagK',meshHierarchy_(1).numDOFs,1);

	%% Assemble&Factorize Stiffness Matrix on Coarsest Level
	[rowIndice, colIndice, ~] = find(ones(24));	
	sK = zeros(24^2, meshHierarchy_(end).numElements);
	for ii=1:meshHierarchy_(end).numElements
		sK(:,ii) = reshape(Ks(:,:,ii), 24^2, 1);
	end
	eNodMat = meshHierarchy_(end).eNodMat;
	eDofMat = [3*eNodMat-2 3*eNodMat-1 3*eNodMat];
	eDofMat = eDofMat(:,reOrdering);
	iK = eDofMat(:,rowIndice);
	jK = eDofMat(:,colIndice);
	KcoarsestLevel = sparse(iK, jK, sK');
	[cholFac_, ~, cholPermut_] = chol(KcoarsestLevel(meshHierarchy_(end).freeDOFs, meshHierarchy_(end).freeDOFs),'lower');
	
	clear Ks KsPrevious
end
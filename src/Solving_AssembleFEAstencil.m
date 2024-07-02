function Solving_AssembleFEAstencil()
	global meshHierarchy_;
	global numLevels_;
	global cholFac_; global cholPermut_;

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
		numProjectDOFs = (spanWidth+1)^3*3;
		meshHierarchy_(ii).storingState = 1;
		meshHierarchy_(ii).Ke = meshHierarchy_(ii-1).Ke*spanWidth;
		numElements = meshHierarchy_(ii).numElements;
		Ks = repmat(meshHierarchy_(ii).Ke, 1,1,numElements);					
		finerKes = zeros(24*24,spanWidth^3);
		elementUpwardMap = meshHierarchy_(ii).elementUpwardMap;			
		if 2==ii
			if size(meshHierarchy_(1).Ks,3)==meshHierarchy_(1).numElements 
				%%for temporary test of anisotropic material law
				KsPrevious = meshHierarchy_(ii-1).Ks;
				for jj=1:numElements
					iFinerEles = elementUpwardMap(jj,:);
					solidEles = find(0~=iFinerEles);
					iFinerEles = iFinerEles(solidEles);
					sK = finerKes;
					tarKes = KsPrevious(:,:,iFinerEles);
					for kk=1:length(solidEles)
						sK(:,solidEles(kk)) = reshape(tarKes(:,:,kk),24^2,1);
					end
					tmpK = sparse(iK, jK, sK, numProjectDOFs, numProjectDOFs);
					tmpK = interpolatingKe' * tmpK * interpolatingKe;
					Ks(:,:,jj) = full(tmpK);
				end								
			else
				iKe = meshHierarchy_(ii-1).Ke;
				iKs = reshape(iKe, 24*24, 1);
				eleModulus = meshHierarchy_(1).eleModulus;
				for jj=1:numElements
					sonEles = elementUpwardMap(jj,:);
					solidEles = find(0~=sonEles);
					sK = finerKes;
					sK(:,solidEles) = iKs .* eleModulus(sonEles(solidEles));
					tmpK = sparse(iK, jK, sK, numProjectDOFs, numProjectDOFs);
					tmpK = interpolatingKe' * tmpK * interpolatingKe;
					Ks(:,:,jj) = full(tmpK);							
				end						
			end	
		else
			KsPrevious = meshHierarchy_(ii-1).Ks;
			for jj=1:numElements
				iFinerEles = elementUpwardMap(jj,:);
				solidEles = find(0~=iFinerEles);
				iFinerEles = iFinerEles(solidEles);
				sK = finerKes;
				tarKes = KsPrevious(:,:,iFinerEles);
				for kk=1:length(solidEles)
					sK(:,solidEles(kk)) = reshape(tarKes(:,:,kk),24^2,1);
				end
				tmpK = sparse(iK, jK, sK, numProjectDOFs, numProjectDOFs);
				tmpK = interpolatingKe' * tmpK * interpolatingKe;
				Ks(:,:,jj) = full(tmpK);
			end									
		end			
		meshHierarchy_(ii).Ks = Ks;
	end		
	
	%% initialize smoother
	for ii=1:length(meshHierarchy_)-1
		diagK = zeros(meshHierarchy_(ii).numDOFs,1);
		% eDofMat = meshHierarchy_(ii).eDofMat;
		eDofMat = [3*meshHierarchy_(ii).eNodMat-2 3*meshHierarchy_(ii).eNodMat-1 3*meshHierarchy_(ii).eNodMat];
		eDofMat = eDofMat(:, reOrdering);
		numElements = meshHierarchy_(ii).numElements;
		Ks = meshHierarchy_(ii).Ks;
		if 1==size(Ks,3)
			diagKe = diag(meshHierarchy_(ii).Ks);
			eleModulus = meshHierarchy_(ii).eleModulus;
			for jj=1:numElements
				dofIndex = eDofMat(jj,:)';
				diagK(dofIndex) = diagK(dofIndex) + diagKe*eleModulus(jj);
			end								
		else
			if 1==ii
				eleModulus = meshHierarchy_(ii).eleModulus;
			else
				eleModulus = ones(1,numElements);
			end
			for jj=1:numElements
				dofIndex = eDofMat(jj,:)';
				Ke = Ks(:,:,jj);
				diagK(dofIndex) = diagK(dofIndex) + diag(Ke)*eleModulus(jj);
			end							
		end		
		meshHierarchy_(ii).diagK = diagK;			
	end	

	%% Assemble&Factorize Stiffness Matrix on Coarsest Level
	[rowIndice, colIndice, ~] = find(ones(24));	
	sK = zeros(24^2, meshHierarchy_(end).numElements);
	for ii=1:meshHierarchy_(end).numElements
		sK(:,ii) = reshape(meshHierarchy_(end).Ks(:,:,ii), 24^2, 1);
	end
	eDofMat = [3*meshHierarchy_(end).eNodMat-2 3*meshHierarchy_(end).eNodMat-1 3*meshHierarchy_(end).eNodMat];
	eDofMat = eDofMat(:,reOrdering);
	iK = eDofMat(:,rowIndice);
	jK = eDofMat(:,colIndice);
	KcoarsestLevel = sparse(iK, jK, sK');
	[cholFac_, ~, cholPermut_] = chol(KcoarsestLevel(meshHierarchy_(end).freeDOFs, meshHierarchy_(end).freeDOFs),'lower');
		
	% iK = meshHierarchy_(end).eDofMat(:,rowIndice);
	% jK = meshHierarchy_(end).eDofMat(:,colIndice);
	% KcoarsestLevel = sparse(iK, jK, sK');
	% numBCs = numel(meshHierarchy_(end).freeDOFs);
	% flexibleArr = struct('arr', []);
	% cholFac_ = repmat(flexibleArr, numBCs, 1);
	% cholPermut_ = repmat(flexibleArr, numBCs, 1);
	% for kk=1:numBCs
		% iK = KcoarsestLevel(meshHierarchy_(end).freeDOFs(kk).arr, meshHierarchy_(end).freeDOFs(kk).arr);
		% [cholFac_(kk).arr, ~, cholPermut_(kk).arr] = chol(iK,'lower');
	% end	
end
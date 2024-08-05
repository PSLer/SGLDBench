function productMV = Solving_KbyU_MatrixFree_previous(uVec, varargin)	
	global meshHierarchy_;
	if 1==nargin, iLevel = 1; else, iLevel = varargin{1}; end
	productMV = zeros(meshHierarchy_(iLevel).numNodes,3);
	Ks = meshHierarchy_(iLevel).Ks;
	impOpt = 0; %% 1==previous
	uVec = reshape(uVec,3,meshHierarchy_(iLevel).numNodes)';
	if 1==size(Ks,3)
		if 1~=iLevel, error('Wrong FEA Computing Stencil!'); end
		blockIndex = Solving_MissionPartition(meshHierarchy_(iLevel).numElements, 1.0e7);		
		for jj=1:size(blockIndex,1)	
			rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))'; %%To avoid super-large data block
			iElesNodMat = meshHierarchy_(iLevel).eNodMatHalf(rangeIndex,:);
			iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
			iIntermediateModulus = meshHierarchy_(iLevel).eleModulus(1,rangeIndex);
			subDisVec = zeros(size(iElesNodMat,1),24);
			if impOpt
				tmp = uVec(:,1); subDisVec(:,1:3:24) = tmp(iElesNodMat);
				tmp = uVec(:,2); subDisVec(:,2:3:24) = tmp(iElesNodMat);
				tmp = uVec(:,3); subDisVec(:,3:3:24) = tmp(iElesNodMat);			
				subDisVec = subDisVec*Ks .* iIntermediateModulus(:);
				tmp = subDisVec(:,1:3:24);
				productMV(:,1) = productMV(:,1) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
				tmp = subDisVec(:,2:3:24);
				productMV(:,2) = productMV(:,2) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
				tmp = subDisVec(:,3:3:24);
				productMV(:,3) = productMV(:,3) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);			
			else				
				tmp = uVec(:,1); subDisVec(:,1:3:24) = Vector2Matrix_Indexing_mex(tmp, iElesNodMat); 
				tmp = uVec(:,2); subDisVec(:,2:3:24) = Vector2Matrix_Indexing_mex(tmp, iElesNodMat);
				tmp = uVec(:,3); subDisVec(:,3:3:24) = Vector2Matrix_Indexing_mex(tmp, iElesNodMat);	
				subDisVec = subDisVec*Ks .* iIntermediateModulus(:);
				iElesNodMat = iElesNodMat(:);
				tmp = subDisVec(:,1:3:24);
				productMV(:,1) = productMV(:,1) + Accumarray_mex(iElesNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
				tmp = subDisVec(:,2:3:24);
				productMV(:,2) = productMV(:,2) + Accumarray_mex(iElesNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
				tmp = subDisVec(:,3:3:24);
				productMV(:,3) = productMV(:,3) + Accumarray_mex(iElesNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);		
			end	
		end
	else
		if 1==iLevel
			eleModulus = meshHierarchy_(1).eleModulus;
		else
			eleModulus = ones(1,meshHierarchy_(iLevel).numElements);
		end
		subDisVec = zeros(meshHierarchy_(iLevel).numElements,24);
		eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(iLevel).eNodMatHalf);
		if impOpt
			tmp = uVec(:,1); subDisVec(:,1:3:24) = tmp(eNodMat);
			tmp = uVec(:,2); subDisVec(:,2:3:24) = tmp(eNodMat);
			tmp = uVec(:,3); subDisVec(:,3:3:24) = tmp(eNodMat);		
			for ii=1:meshHierarchy_(iLevel).numElements
				subDisVec(ii,:) = subDisVec(ii,:)*meshHierarchy_(iLevel).Ks(:,:,ii)*eleModulus(ii);
			end
			tmp = subDisVec(:,1:3:24);
			productMV(:,1) = productMV(:,1) + accumarray(eNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,2:3:24);
			productMV(:,2) = productMV(:,2) + accumarray(eNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,3:3:24);
			productMV(:,3) = productMV(:,3) + accumarray(eNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);		
		else
			tmp = uVec(:,1); subDisVec(:,1:3:24) = Vector2Matrix_Indexing_mex(tmp, eNodMat); 
			tmp = uVec(:,2); subDisVec(:,2:3:24) = Vector2Matrix_Indexing_mex(tmp, eNodMat);
			tmp = uVec(:,3); subDisVec(:,3:3:24) = Vector2Matrix_Indexing_mex(tmp, eNodMat);
			for ii=1:meshHierarchy_(iLevel).numElements
				subDisVec(ii,:) = subDisVec(ii,:)*meshHierarchy_(iLevel).Ks(:,:,ii)*eleModulus(ii);
			end
			eNodMat = eNodMat(:);
			tmp = subDisVec(:,1:3:24);
			productMV(:,1) = productMV(:,1) + Accumarray_mex(eNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,2:3:24);
			productMV(:,2) = productMV(:,2) + Accumarray_mex(eNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,3:3:24);
			productMV(:,3) = productMV(:,3) + Accumarray_mex(eNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
		end
	end
	productMV = productMV'; productMV = productMV(:);
	productMV(meshHierarchy_(iLevel).fixedDOFs,1) = 0;	%%Change to Logical Indexing	
end

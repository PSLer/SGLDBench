function productMV = Solving_KbyU_MatrixFree(uVec, varargin)	
	global meshHierarchy_;
	if 1==nargin, iLevel = 1; else, iLevel = varargin{1}; end
	productMV = zeros(meshHierarchy_(iLevel).numNodes,3);
	Ks = meshHierarchy_(iLevel).Ks;
	impOpt = 0; %%1==Previous, 0==New
	if 1==size(Ks,3)
		if 1~=iLevel, error('Wrong FEA Computing Stencil!'); end
		blockIndex = Solving_MissionPartition(meshHierarchy_(iLevel).numElements, 1.0e7);		
		for jj=1:size(blockIndex,1)	
			rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))';
			iElesNodMat = meshHierarchy_(iLevel).eNodMat(rangeIndex,:);
			iIntermediateModulus = meshHierarchy_(iLevel).eleModulus(1,rangeIndex);
			subDisVec = zeros(size(iElesNodMat,1),24);
			tmp = uVec(1:3:end,1); subDisVec(:,1:3:24) = tmp(iElesNodMat);
			tmp = uVec(2:3:end,1); subDisVec(:,2:3:24) = tmp(iElesNodMat);
			tmp = uVec(3:3:end,1); subDisVec(:,3:3:24) = tmp(iElesNodMat);
			subDisVec = subDisVec*Ks .* iIntermediateModulus(:);
			if impOpt %%Previous
				for ii=1:8
					tmp = iElesNodMat(:,ii);
					productMV(tmp,1) = productMV(tmp,1) + subDisVec(:,3*(ii-1)+1);
					productMV(tmp,2) = productMV(tmp,2) + subDisVec(:,3*(ii-1)+2);
					productMV(tmp,3) = productMV(tmp,3) + subDisVec(:,3*(ii-1)+3);
				end				
			else
				tmp = subDisVec(:,1:3:24);
				productMV(:,1) = productMV(:,1) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
				tmp = subDisVec(:,2:3:24);
				productMV(:,2) = productMV(:,2) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
				tmp = subDisVec(:,3:3:24);
				productMV(:,3) = productMV(:,3) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);				
			end		
		end
	else
		if 1==iLevel
			eleModulus = meshHierarchy_(1).eleModulus;
		else
			eleModulus = ones(1,meshHierarchy_(iLevel).numElements);
		end
		
		subDisVec = zeros(meshHierarchy_(iLevel).numElements,24);
		tmp = uVec(1:3:end,1); subDisVec(:,1:3:24) = tmp(meshHierarchy_(iLevel).eNodMat);
		tmp = uVec(2:3:end,1); subDisVec(:,2:3:24) = tmp(meshHierarchy_(iLevel).eNodMat);
		tmp = uVec(3:3:end,1); subDisVec(:,3:3:24) = tmp(meshHierarchy_(iLevel).eNodMat);		
		for ii=1:meshHierarchy_(iLevel).numElements
			subDisVec(ii,:) = subDisVec(ii,:)*meshHierarchy_(iLevel).Ks(:,:,ii)*eleModulus(ii);
		end
		if impOpt
			for ii=1:8
				tmp = meshHierarchy_(iLevel).eNodMat(:,ii);
				productMV(tmp,1) = productMV(tmp,1) + subDisVec(:,3*(ii-1)+1);
				productMV(tmp,2) = productMV(tmp,2) + subDisVec(:,3*(ii-1)+2);
				productMV(tmp,3) = productMV(tmp,3) + subDisVec(:,3*(ii-1)+3);
			end			
		else
			tmp = subDisVec(:,1:3:24);
			productMV(:,1) = productMV(:,1) + accumarray(meshHierarchy_(iLevel).eNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,2:3:24);
			productMV(:,2) = productMV(:,2) + accumarray(meshHierarchy_(iLevel).eNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,3:3:24);
			productMV(:,3) = productMV(:,3) + accumarray(meshHierarchy_(iLevel).eNodMat(:),tmp(:),[meshHierarchy_(iLevel).numNodes 1]);		
		end
	
	end
	productMV = productMV'; productMV = productMV(:);
	productMV(meshHierarchy_(iLevel).fixedDOFs,1) = 0;		
end

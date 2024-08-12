function productMV = Solving_KbyU_MatrixFree(uVec)	
	global meshHierarchy_;
	global MEXfunc_;
	iLevel = 1;
	blockSize = 1.0e7;
	Ks = meshHierarchy_(1).Ks;
	eleModulus = meshHierarchy_(1).eleModulus;
	if MEXfunc_
		productMV = Solving_KbyU_MatrixFree_mex(uVec, meshHierarchy_(1).eNodMat, Ks, eleModulus(:), blockSize);
	else %% full mex
		uVec = reshape(uVec,3,meshHierarchy_(1).numNodes)';
		productMV = zeros(meshHierarchy_(1).numNodes,3);
		%%To avoid super-large data block
        blockIndex = Solving_MissionPartition(meshHierarchy_(iLevel).numElements, blockSize);
		for jj=1:size(blockIndex,1)
		
			if 1==size(blockIndex,1)
				% iElesNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(1).eNodMatHalf);
				iElesNodMat = meshHierarchy_(1).eNodMat;
				iIntermediateModulus = eleModulus;
			else
				rangeIndex = (blockIndex(jj,1):blockIndex(jj,2));
				% iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
				% iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
				iElesNodMat = meshHierarchy_(1).eNodMat(rangeIndex,:);
				iIntermediateModulus = eleModulus(1,rangeIndex);
			end
			subDisVec = zeros(size(iElesNodMat,1),24);
			if MEXfunc_		
				tmp = uVec(:,1); subDisVec(:,1:3:24) = Vector2Matrix_Indexing_mex_openMP(tmp, iElesNodMat); 
				tmp = uVec(:,2); subDisVec(:,2:3:24) = Vector2Matrix_Indexing_mex_openMP(tmp, iElesNodMat);
				tmp = uVec(:,3); subDisVec(:,3:3:24) = Vector2Matrix_Indexing_mex_openMP(tmp, iElesNodMat);	
				subDisVec = subDisVec*Ks .* iIntermediateModulus(:);					
				iElesNodMat = iElesNodMat(:);
				tmp = subDisVec(:,1:3:24);
				productMV(:,1) = productMV(:,1) + Accumarray_mex_openMP(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);
				tmp = subDisVec(:,2:3:24);
				productMV(:,2) = productMV(:,2) + Accumarray_mex_openMP(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);
				tmp = subDisVec(:,3:3:24);
				productMV(:,3) = productMV(:,3) + Accumarray_mex_openMP(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);			
			else				
				tmp = uVec(:,1); subDisVec(:,1:3:24) = tmp(iElesNodMat);
				tmp = uVec(:,2); subDisVec(:,2:3:24) = tmp(iElesNodMat);
				tmp = uVec(:,3); subDisVec(:,3:3:24) = tmp(iElesNodMat);			
				subDisVec = subDisVec*Ks .* iIntermediateModulus(:);
				tmp = subDisVec(:,1:3:24);
				productMV(:,1) = productMV(:,1) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(1).numNodes 1]);
				tmp = subDisVec(:,2:3:24);
				productMV(:,2) = productMV(:,2) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(1).numNodes 1]);
				tmp = subDisVec(:,3:3:24);
				productMV(:,3) = productMV(:,3) + accumarray(iElesNodMat(:),tmp(:),[meshHierarchy_(1).numNodes 1]);					
			end	
		end
		productMV = productMV'; productMV = productMV(:);		
	end
	productMV(meshHierarchy_(1).fixedDOFs,1) = 0;	%%Change to Logical Indexing
end

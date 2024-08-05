function productMV = Solving_KbyU_MatrixFree(uVec, varargin)	
	global meshHierarchy_;
	global missionPartitionIndexing_;
	productMV = zeros(meshHierarchy_(1).numNodes,3);
	Ks = meshHierarchy_(1).Ks;
	impOpt = 0; %% 1==previous
	uVec = reshape(uVec,3,meshHierarchy_(1).numNodes)';
	%%To avoid super-large data block
	for jj=1:numel(missionPartitionIndexing_,1)			
		rangeIndex = missionPartitionIndexing_(jj).logicalIndexingElement;
		iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
		iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
		iIntermediateModulus = meshHierarchy_(1).eleModulus(1,rangeIndex);
		subDisVec = zeros(size(iElesNodMat,1),24);
		if impOpt
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
		else				
			tmp = uVec(:,1); subDisVec(:,1:3:24) = Vector2Matrix_Indexing_mex(tmp, iElesNodMat); 
			tmp = uVec(:,2); subDisVec(:,2:3:24) = Vector2Matrix_Indexing_mex(tmp, iElesNodMat);
			tmp = uVec(:,3); subDisVec(:,3:3:24) = Vector2Matrix_Indexing_mex(tmp, iElesNodMat);	
			subDisVec = subDisVec*Ks .* iIntermediateModulus(:);
			iElesNodMat = iElesNodMat(:);
			tmp = subDisVec(:,1:3:24);
			productMV(:,1) = productMV(:,1) + Accumarray_mex(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);
			tmp = subDisVec(:,2:3:24);
			productMV(:,2) = productMV(:,2) + Accumarray_mex(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);
			tmp = subDisVec(:,3:3:24);
			productMV(:,3) = productMV(:,3) + Accumarray_mex(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);		
		end	
	end
	productMV = productMV'; productMV = productMV(:);
	productMV(meshHierarchy_(1).fixedDOFs,1) = 0;	%%Change to Logical Indexing	
end

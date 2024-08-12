function productMV = Solving_KbyU_MatrixFree_milestone(uVec, varargin)	
	global meshHierarchy_;
	global MEXfunc_;
	if 1==nargin, iLevel = 1; else, iLevel = varargin{1}; end
% tStart1 = tic;	
	productMV = zeros(meshHierarchy_(1).numNodes,3);
	Ks = meshHierarchy_(1).Ks;
	uVec = reshape(uVec,3,meshHierarchy_(1).numNodes)';
% tStart1Total = toc(tStart1);	
% tStart2Total = 0;
% tStart3Total = 0;
% tStart4Total = 0;
% tStart5Total = 0;
% tStart6Total = 0;
	if 1==size(Ks,3)
		%%To avoid super-large data block
        blockIndex = Solving_MissionPartition(meshHierarchy_(iLevel).numElements, 1.0e7);
		for jj=1:size(blockIndex,1)
% tStart2 = tic;			
			if 1==size(blockIndex,1)
				% iElesNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(1).eNodMatHalf);
				iElesNodMat = meshHierarchy_(1).eNodMat;
				iIntermediateModulus = meshHierarchy_(1).eleModulus;
			else
				rangeIndex = (blockIndex(jj,1):blockIndex(jj,2));
				% iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
				% iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
				iElesNodMat = meshHierarchy_(1).eNodMat(rangeIndex,:);
				iIntermediateModulus = meshHierarchy_(1).eleModulus(1,rangeIndex);
			end
% tStart2Total = tStart2Total + toc(tStart2);			
% tStart3 = tic;
			subDisVec = zeros(size(iElesNodMat,1),24);
			if MEXfunc_		
				tmp = uVec(:,1); subDisVec(:,1:3:24) = Vector2Matrix_Indexing_mex_openMP(tmp, iElesNodMat); 
				tmp = uVec(:,2); subDisVec(:,2:3:24) = Vector2Matrix_Indexing_mex_openMP(tmp, iElesNodMat);
				tmp = uVec(:,3); subDisVec(:,3:3:24) = Vector2Matrix_Indexing_mex_openMP(tmp, iElesNodMat);
% tStart3Total = tStart3Total + toc(tStart3);
% tStart4 = tic;	
				subDisVec = subDisVec*Ks .* iIntermediateModulus(:);
% tStart4Total = tStart4Total + toc(tStart4);
% tStart5 = tic;					
				iElesNodMat = iElesNodMat(:);
				tmp = subDisVec(:,1:3:24);
				productMV(:,1) = productMV(:,1) + Accumarray_mex_openMP(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);
				tmp = subDisVec(:,2:3:24);
				productMV(:,2) = productMV(:,2) + Accumarray_mex_openMP(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);
				tmp = subDisVec(:,3:3:24);
				productMV(:,3) = productMV(:,3) + Accumarray_mex_openMP(iElesNodMat,tmp(:),[meshHierarchy_(1).numNodes 1]);
% tStart5Total = tStart5Total + toc(tStart5);				
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
	else
		if 1==iLevel
			eleModulus = meshHierarchy_(1).eleModulus;
		else
			eleModulus = ones(1,meshHierarchy_(iLevel).numElements);
		end
		subDisVec = zeros(meshHierarchy_(iLevel).numElements,24);
		% eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(iLevel).eNodMatHalf);
		eNodMat = meshHierarchy_(iLevel).eNodMat;
		if MEXfunc_
			tmp = uVec(:,1); subDisVec(:,1:3:24) = Vector2Matrix_Indexing_mex(tmp, eNodMat); 
			tmp = uVec(:,2); subDisVec(:,2:3:24) = Vector2Matrix_Indexing_mex(tmp, eNodMat);
			tmp = uVec(:,3); subDisVec(:,3:3:24) = Vector2Matrix_Indexing_mex(tmp, eNodMat);
			for ii=1:meshHierarchy_(iLevel).numElements
				subDisVec(ii,:) = subDisVec(ii,:)*meshHierarchy_(iLevel).Ks(:,:,ii);
			end
			eNodMat = eNodMat(:);
			tmp = subDisVec(:,1:3:24);
			productMV(:,1) = productMV(:,1) + Accumarray_mex(eNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,2:3:24);
			productMV(:,2) = productMV(:,2) + Accumarray_mex(eNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);
			tmp = subDisVec(:,3:3:24);
			productMV(:,3) = productMV(:,3) + Accumarray_mex(eNodMat,tmp(:),[meshHierarchy_(iLevel).numNodes 1]);		
		else		
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
		end		
	end
% tStart6 = tic;	
	productMV = productMV'; productMV = productMV(:);
	productMV(meshHierarchy_(1).fixedDOFs,1) = 0;	%%Change to Logical Indexing
% tStart6Total = tStart6Total + toc(tStart6);
% disp(['tStart1Total: Initialization Costs: ', sprintf('%10.3g',tStart1Total) 's']);
% disp(['tStart2Total: Recover eNodMat Costs: ', sprintf('%10.3g',tStart2Total) 's']);
% disp(['tStart3Total: Vector-to-Matrix Indexing Costs: ', sprintf('%10.3g',tStart3Total) 's']);
% disp(['tStart4Total: Square Matrix Times Band Matrix Costs: ', sprintf('%10.3g',tStart4Total) 's']);
% disp(['tStart5Total: Matrix-to-Vector Indexing Costs: ', sprintf('%10.3g',tStart5Total) 's']);
% disp(['tStart6Total: Finilize Output Costs: ', sprintf('%10.3g',tStart6Total) 's']);
% disp(['Total: ', sprintf('%10.3g',tStart1Total + tStart2Total + tStart3Total + tStart4Total + tStart5Total + tStart6Total) 's']);	
end

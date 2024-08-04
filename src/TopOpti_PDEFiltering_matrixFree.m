function tar = TopOpti_PDEFiltering_matrixFree(src)
	global meshHierarchy_;
	global diagKePDE_;
	%%Element to Node
	src = double(src);
	tmpVal = zeros(meshHierarchy_(1).numNodes,1);
if 1
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,1)+1) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,2)+1) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,3)+1) + src*(1/8);
	tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)+1) = tmpVal(meshHierarchy_(1).eNodMatHalf(:,4)+1) + src*(1/8);
else
	
end
	% src = double(tmpVal);	
	src = tmpVal;
	
	%% Solving on Node
	PtV = @(x) diagKePDE_ .* x;
	%tar = pcg(@MatTimesVec_matrixFree, src, 1.0e-6, 200, PtV);
	tar = Solving_PreconditionedConjugateGradientSolver(@MatTimesVec_matrixFree, PtV, src, 1.0e-6, 200, 'printP_ON');
	
	%%Node to Element
	% tar = single(tar);
	tmpVal = zeros(meshHierarchy_(1).numElements,1);
	blockIndex = Solving_MissionPartition(meshHierarchy_(1).numElements, 1.0e7);
	for jj=1:size(blockIndex,1)	
		rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))'; %%To avoid super-large data block
		iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
		iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
		tmpVal(rangeIndex,1) = sum(tar(iElesNodMat),2);
	end
	tar = tmpVal*(1/8);	
end


function productMV = MatTimesVec_matrixFree(uVec)	
	global meshHierarchy_;
	global KePDE_;
	productMV = zeros(meshHierarchy_(1).numNodes,1);
	Ks = KePDE_;
if 1
	blockIndex = Solving_MissionPartition(meshHierarchy_(1).numElements, 1.0e7);		
	for jj=1:size(blockIndex,1)	
		rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))'; %%To avoid super-large data block
		iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
		iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
		subDisVec = uVec(iElesNodMat);
		subDisVec = subDisVec*Ks;
		productMV = productMV + accumarray(iElesNodMat(:),subDisVec(:),[meshHierarchy_(1).numNodes 1]);
	end
else
	iElesNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(1).eNodMatHalf)';
	subDisVec = uVec(iElesNodMat);
	parpool('Threads');
	iKs = Ks(1,:); iNodes = iElesNodMat(:,1);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii); productMV(iNodes(ii)) = productMV(iNodes(ii)) + iVal;
	% end
	% iKs = Ks(2,:); iNodes = iElesNodMat(:,2);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii);
		% tmpNode = iNodes(ii);
		% productMV(tmpNode) = iVal;
	% end
	% iKs = Ks(3,:); iNodes = iElesNodMat(:,3);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii);
		% tmpNode = iNodes(ii);
		% iVal2 = productMV(tmpNode);
		% productMV(tmpNode) = iVal2 + iVal;
	% end	
	% iKs = Ks(4,:); iNodes = iElesNodMat(:,4);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii);
		% tmpNode = iNodes(ii);
		% iVal2 = productMV(tmpNode);
		% productMV(tmpNode) = iVal2 + iVal;
	% end
	% iKs = Ks(5,:); iNodes = iElesNodMat(:,5);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii);
		% tmpNode = iNodes(ii);
		% iVal2 = productMV(tmpNode);
		% productMV(tmpNode) = iVal2 + iVal;
	% end	
	% iKs = Ks(6,:); iNodes = iElesNodMat(:,6);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii);
		% tmpNode = iNodes(ii);
		% iVal2 = productMV(tmpNode);
		% productMV(tmpNode) = iVal2 + iVal;
	% end	
	% iKs = Ks(7,:); iNodes = iElesNodMat(:,7);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii);
		% tmpNode = iNodes(ii);
		% iVal2 = productMV(tmpNode);
		% productMV(tmpNode) = iVal2 + iVal;
	% end	
	% iKs = Ks(8,:); iNodes = iElesNodMat(:,8);
	% parfor ii=1:meshHierarchy_(1).numElements
		% iVal = iKs * subDisVec(:,ii);
		% tmpNode = iNodes(ii);
		% iVal2 = productMV(tmpNode);
		% productMV(tmpNode) = iVal2 + iVal;
	% end	
	% delete(gcp('nocreate'));
end
end

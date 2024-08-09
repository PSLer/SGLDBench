function tar = TopOpti_PDEFiltering_matrixFree(src)
	global meshHierarchy_;
	global diagKePDE_;
	%%Element to Node
	src = double(src);
	tmpVal = zeros(meshHierarchy_(1).numNodes,1);
	values = src(:)*(1/8);
	
	%%Ref
	for jj=1:8
		tmpVal = tmpVal + accumarray(meshHierarchy_(1).eNodMat(:,jj), values, [meshHierarchy_(1).numNodes, 1]);
		%%Ref
		% tmpVal(meshHierarchy_(1).eNodMat(:,jj),1) = tmpVal(meshHierarchy_(1).eNodMat(:,jj),1) + values;
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
	blockIndex = Solving_MissionPartition(meshHierarchy_(1).numElements, 3.0e7);
	for jj=1:size(blockIndex,1)	
		rangeIndex = (blockIndex(jj,1):blockIndex(jj,2))'; %%To avoid super-large data block
		iElesNodMat = meshHierarchy_(1).eNodMat(rangeIndex,:);
		tmpVal(rangeIndex,1) = sum(tar(iElesNodMat),2);
	end
	tar = tmpVal*(1/8);
end

function productMV = MatTimesVec_matrixFree(uVec)	
	global meshHierarchy_;
	global KePDE_;
	global MEXfunc_;
	blockSize = 3.0e7;
	productMV = zeros(meshHierarchy_(1).numNodes,1);
	Ks = KePDE_;
	if MEXfunc_
		productMV = Solving_KbyU_MatrixFree8x8_mex(uVec, meshHierarchy_(1).eNodMat, Ks, blockSize);
	else
		blockIndex = Solving_MissionPartition(meshHierarchy_(1).numElements, blockSize);
		if 0 %% partly MEX
			for jj=1:size(blockIndex,1)
				iElesNodMat = meshHierarchy_(1).eNodMat((blockIndex(jj,1):blockIndex(jj,2)),:);
				subDisVec = Vector2Matrix_Indexing_mex(uVec, iElesNodMat);
				subDisVec = subDisVec*Ks;
				productMV = productMV + Accumarray_mex(iElesNodMat(:),subDisVec(:),[meshHierarchy_(1).numNodes 1]);
			end		
		else
			for jj=1:size(blockIndex,1)
				iElesNodMat = meshHierarchy_(1).eNodMat((blockIndex(jj,1):blockIndex(jj,2))',:);
				subDisVec = uVec(iElesNodMat);
				subDisVec = subDisVec*Ks;
				productMV = productMV + accumarray(iElesNodMat(:),subDisVec(:),[meshHierarchy_(1).numNodes 1]);
			end		
		end
	end
end

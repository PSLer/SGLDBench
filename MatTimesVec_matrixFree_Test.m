% clc;

global meshHierarchy_;
global KePDE_;
uVec = rand(meshHierarchy_(1).numNodes,1);
eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(1).eNodMatHalf);
blockIndex = [1 meshHierarchy_(1).numElements];
if 0	
	rangeIndex = (blockIndex(1):blockIndex(2))'; %%To avoid super-large data block
else
	rangeIndex = false(meshHierarchy_(1).numElements,1); rangeIndex(blockIndex(1):blockIndex(2),1) = true;
end


tStart0 = tic;
productMV = zeros(meshHierarchy_(1).numNodes,1);
Ks = KePDE_;
tStart1 = tic;
if 0
	iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
	iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
else
	iElesNodMat = eNodMat(rangeIndex,:); iElesNodMat_flat = iElesNodMat(:);
end
disp(['tStart1: Recover eNodMat Costs: ', sprintf('%10.3g',toc(tStart1)) 's']);

tStart2 = tic;
opt = 5; %%1: original; 2: Vectorized; 3: Prelocation Loop; 4: parfor-loop; 5: MEX (single thread); 6: MEX (openMP)
switch opt
	case 1
		subDisVec = uVec(iElesNodMat);
	case 2
		subDisVec = uVec(iElesNodMat_flat);
		subDisVec = reshape(subDisVec, numel(iElesNodMat_flat)/8, 8);	
	case 3
		for ii=1:rows
			for jj=1:cols
				subDisVec(ii,jj) = uVec(iElesNodMat(ii,jj));
			end
		end
	case 4
		[rows, cols] = size(iElesNodMat);
		subDisVec = zeros(rows*cols,1);		
		parpool('Threads');
		parfor ii=1:rows*cols
			subDisVec(ii) = uVec(iElesNodMat_flat(ii));
		end
		delete(gcp('nocreate'));
		subDisVec = reshape(subDisVec, rows, cols);
	case 5
		subDisVec = Vector2Matrix_Indexing_mex(uVec, iElesNodMat);         
	case 6
		subDisVec = Vector2Matrix_Indexing_mex_openMP(uVec, iElesNodMat); 
end
disp(['tStart2: Reshape Vector to Matrix Costs: ', sprintf('%10.3g',toc(tStart2)) 's']);

tStart3 = tic;
subDisVec = subDisVec*Ks;
disp(['tStart3: 8-by-8 Square Matrix Times 8-by-10000000 Stripe Matrix Costs: ', sprintf('%10.3g',toc(tStart3)) 's']);

tStart4 = tic;
if 0
	productMV = productMV + accumarray(iElesNodMat(:),subDisVec(:),[meshHierarchy_(1).numNodes 1]);
else
	productMV = productMV + Accumarray_mex(iElesNodMat(:), subDisVec(:), [meshHierarchy_(1).numNodes, 1]);
end
disp(['tStart4: Project Result in Vector to Matrix Form Costs: ', sprintf('%10.3g',toc(tStart4)) 's']);

disp(['tStart0: Conduct Matrix times Vector in Matrix-free Format Costs (Whole Process): ', sprintf('%10.3g',toc(tStart0)) 's']);

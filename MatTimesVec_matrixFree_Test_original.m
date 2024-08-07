clc;

global meshHierarchy_;
global KePDE_;
uVec = rand(meshHierarchy_(1).numNodes,1);

tStart0 = tic;
productMV = zeros(meshHierarchy_(1).numNodes,1);
Ks = KePDE_;
blockIndex = [1 meshHierarchy_(1).numElements];		
rangeIndex = (blockIndex(1):blockIndex(2))'; %%To avoid super-large data block
iElesNodMat = meshHierarchy_(1).eNodMatHalf(rangeIndex,:);
iElesNodMat = Common_RecoverHalfeNodMat(iElesNodMat);
subDisVec = uVec(iElesNodMat);
subDisVec = subDisVec*Ks;
productMV = productMV + accumarray(iElesNodMat(:),subDisVec(:),[meshHierarchy_(1).numNodes 1]);
disp(['Conduct Matrix times Vector in Matrix-free Format Costs: ', sprintf('%10.3g',toc(tStart0)) 's']);
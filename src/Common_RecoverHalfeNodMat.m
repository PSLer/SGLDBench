function eNodMat = Common_RecoverHalfeNodMat(eNodMatHalf)
	if 4~=size(eNodMatHalf,2), eNodMat = []; return; end
	numEles = size(eNodMatHalf,1);
	eNodMat = zeros(numEles,8,'int32');
	eNodMat(:,[3 4 7 8]) = eNodMatHalf;
	eNodMat(:,[2 1 6 5]) = eNodMatHalf + 1;
end
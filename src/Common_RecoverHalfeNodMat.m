function eNodMat = Common_RecoverHalfeNodMat(eNodMatHalf)
	if 4~=size(eNodMatHalf,2), eNodMat = []; return; end
	numEles = size(eNodMatHalf,1);
	eNodMat = zeros(numEles,8,'int32');
	if 1
		eNodMat(:,[3 4 7 8]) = eNodMatHalf;
		eNodMat(:,[2 1 6 5]) = eNodMatHalf + 1;	
	else
		col1 = false(1,8); col1(1,[3 4 7 8]) = true;
		col2 = false(1,8); col2(1,[2 1 6 5]) = true;
		eNodMat(:,col1) = eNodMatHalf;
		eNodMat(:,col2) = eNodMatHalf + 1;
		eNodMat(:, [1 2 5 6]) = eNodMat(:,[2 1 6 5]);
	end
	
end
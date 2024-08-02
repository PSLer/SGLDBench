function TopOpti_BuildPDEfilter()
	global meshHierarchy_;
	global rHatMin_; 
	global KF_; global TF_; global LF_; global Permut_;
    global PDEfilterSolver_;
    
	rHatMax_ = rHatMin_;
	numElements = meshHierarchy_(1).numElements;
	numNodes = meshHierarchy_(1).numNodes;
	iNumNode = 8;
	eNodMat = Common_RecoverHalfeNodMat(meshHierarchy_(1).eNodMatHalf);
	iTF = reshape(eNodMat,iNumNode*numElements,1);
	jTF = reshape(repmat(1:int32(numElements),iNumNode,1)',iNumNode*numElements,1);
	sTF = repmat(1/iNumNode,iNumNode*numElements,1);
	TF_ = sparse(iTF,jTF,sTF);
	
	iKF = reshape(kron(eNodMat,ones(iNumNode,1,'int32'))',(iNumNode)^2*numElements,1);
	jKF = reshape(kron(eNodMat,ones(1,iNumNode,'int32'))',(iNumNode)^2*numElements,1);
	
	[s, t, p, w] = FEA_GaussianIntegral();
	N = FEA_ShapeFunction(s, t, p);
	dShape = FEA_DeShapeFunction(s,t,p);
	
	%%Jacobian Matrix, corresponding to the commonly used 2x2x2 cubic element in natural coordinate system
	% detJ = ones(8,1);
	% invJ = eye(24,24)*2;
	%%CellSize
	CellSize = 1; %%Always treated as a unit cell
	detJ = CellSize^3 /8 * ones(8,1); %%Sub-Volume

	wgt = w.*detJ;
	KEF0 = dShape'*dShape;
	KEF1 = N'*diag(wgt)*N;

	eleMapForward = meshHierarchy_(1).eleMapForward;
	nelx = meshHierarchy_(1).resX;
	nely = meshHierarchy_(1).resY;
	nelz = meshHierarchy_(1).resZ;
    sKF = KEF1(:)*ones(1,numElements);
	stepSize = (rHatMax_-rHatMin_) * meshHierarchy_(1).eleSize(1) / (nelx-1);
	for ii=1:nelx
		iRmin = (rHatMin_ * meshHierarchy_(1).eleSize(1) + (ii-1)*stepSize)/2/sqrt(3);
		iKEF = iRmin^2*KEF0 + KEF1; iKEF = iKEF(:);
		for jj=1:nely
			for kk=1:nelz
				iEle = eleMapForward((kk-1)*nelx*nely+(ii-1)*nely+jj);
				if iEle
					sKF(:,iEle) = iKEF;
				end							
			end
		end
	end
	sKF = sKF(:);

	blockIndex = Solving_MissionPartition(size(iKF,1), 1.0e8);
	KF_ = sparse(numNodes, numNodes);
	for ii=1:size(blockIndex,1)
		rangeIndex = (blockIndex(ii,1):blockIndex(ii,2))';
		tmpKF = sparse(iKF(rangeIndex,:), jKF(rangeIndex,:), sKF(rangeIndex,:), numNodes, numNodes);					
		KF_ = KF_ + tmpKF;
	end		
	PDEfilterSolver_ = 0; %%Iterative
	LF_ = ichol(KF_);
	clearvars iKF jKF sKF
end
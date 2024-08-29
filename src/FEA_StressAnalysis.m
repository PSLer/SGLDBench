function [cartesianStressField, vonMisesStressField] = FEA_StressAnalysis()
	%% sigma_xx, sigma_yy, sigma_zz, tadisyz, tadiszx, tadisxy (3D)	
	global meshHierarchy_;	
	global U_; 
	global cellSize_;

	cartesianStressField = [];
	vonMisesStressField = [];
	
	if isempty(U_), return; end

	%%Element Strain Matrix (eleB)
	%%Jacobian Matrix, corresponding to the commonly used 2x2x2 cubic element in natural coordinate system
	% detJ = ones(8,1);
	% invJ = eye(24,24)*2;
	%%CellSize	
	detJ = cellSize_^3 /8 * ones(8,1); %%Sub-Volume
	invJ = eye(24,24)*2*(1/cellSize_);
	[s, t, p, ~] = FEA_GaussianIntegral();
	dShape = FEA_DeShapeFunction(s,t,p);
	eleB = FEA_ElementStrainMatrix(dShape, invJ);

	%%Element Elasticity Matrix (D-mat)
	S = FEA_HookeLaw();
	eleD = FEA_ElementElasticityMatrix(S);
	U_ = reshape(U_, 3, meshHierarchy_(1).numNodes)';
	
	%%Cartesian Stress
    cartesianStressField = zeros(meshHierarchy_(1).numNodes, 6);
	OTP = OuterInterpolationMat();
	eleModulus = meshHierarchy_(1).eleModulus;
	% eNodMatHalf = meshHierarchy_(1).eNodMatHalf;
	eNodMat = meshHierarchy_(1).eNodMat;
	for ii=1:meshHierarchy_(1).numElements
		% relativeNodesIndex = eNodMatHalf(ii,:);
		% relativeNodesIndex = Common_RecoverHalfeNodMat(relativeNodesIndex);
		relativeNodesIndex = eNodMat(ii,:);
		iEleU = U_(relativeNodesIndex,:)'; iEleU = iEleU(:);
		cartesianStressOnGaussIntegralPoints = eleModulus(ii)*eleD * (eleB*iEleU);	
		midVar = OTP*cartesianStressOnGaussIntegralPoints;
		midVar = reshape(midVar, 6, 8)';	
		cartesianStressField(relativeNodesIndex,:) = midVar + cartesianStressField(relativeNodesIndex,:);
	end
	cartesianStressField = cartesianStressField./meshHierarchy_(1).numNod2ElesVec;
	U_ = reshape(U_', numel(U_), 1);
	
	%%Von Mises Stress
	vonMisesStressField = sqrt(0.5*((cartesianStressField(:,1)-cartesianStressField(:,2)).^2 + ...
		(cartesianStressField(:,2)-cartesianStressField(:,3)).^2 + (cartesianStressField(:,3)...
			-cartesianStressField(:,1)).^2 ) + 3*( cartesianStressField(:,6).^2 + cartesianStressField(:,4).^2 + ...
				cartesianStressField(:,5).^2 ));		
end

function outerInterpolationMatrix = OuterInterpolationMat()
	[s, t, p, ~] = FEA_GaussianIntegral();
	N = FEA_ShapeFunction(s,t,p);
	sFM = sparse(48,48);
	ii = 6*(1:8);
	sFM(1,ii-5) = N(1,:); sFM(2,ii-4) = N(1,:); sFM(3,ii-3) = N(1,:);
	sFM(4,ii-2) = N(1,:); sFM(5,ii-1) = N(1,:); sFM(6,ii) = N(1,:);
	
	sFM(7,ii-5) = N(2,:); sFM(8,ii-4) = N(2,:); sFM(9,ii-3) = N(2,:);
	sFM(10,ii-2) = N(2,:); sFM(11,ii-1) = N(2,:); sFM(12,ii) = N(2,:);

	sFM(13,ii-5) = N(3,:); sFM(14,ii-4) = N(3,:); sFM(15,ii-3) = N(3,:);
	sFM(16,ii-2) = N(3,:); sFM(17,ii-1) = N(3,:); sFM(18,ii) = N(3,:);	

	sFM(19,ii-5) = N(4,:); sFM(20,ii-4) = N(4,:); sFM(21,ii-3) = N(4,:);
	sFM(22,ii-2) = N(4,:); sFM(23,ii-1) = N(4,:); sFM(24,ii) = N(4,:);

	sFM(25,ii-5) = N(5,:); sFM(26,ii-4) = N(5,:); sFM(27,ii-3) = N(5,:);
	sFM(28,ii-2) = N(5,:); sFM(29,ii-1) = N(5,:); sFM(30,ii) = N(5,:);	

	sFM(31,ii-5) = N(6,:); sFM(32,ii-4) = N(6,:); sFM(33,ii-3) = N(6,:);
	sFM(34,ii-2) = N(6,:); sFM(35,ii-1) = N(6,:); sFM(36,ii) = N(6,:);	

	sFM(37,ii-5) = N(7,:); sFM(38,ii-4) = N(7,:); sFM(39,ii-3) = N(7,:);
	sFM(40,ii-2) = N(7,:); sFM(41,ii-1) = N(7,:); sFM(42,ii) = N(7,:);	

	sFM(43,ii-5) = N(8,:); sFM(44,ii-4) = N(8,:); sFM(45,ii-3) = N(8,:);
	sFM(46,ii-2) = N(8,:); sFM(47,ii-1) = N(8,:); sFM(48,ii) = N(8,:);		
	outerInterpolationMatrix = inv(sFM);
end

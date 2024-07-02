function Ke = FEA_VoxelBasedElementStiffnessMatrix()
	global cellSize_;
	%% Compute the Element Stiffness Matrix of the 8-node Cubic Finite Element with 8 Gaussian Integration Points
	
	%%Element Strain Matrix (eleB)
	%%Jacobian Matrix, corresponding to the commonly used 2x2x2 cubic element in natural coordinate system
	% detJ = ones(8,1);
	% invJ = eye(24,24)*2;
	%%CellSize
	detJ = cellSize_^3 /8 * ones(8,1); %%Sub-Volume
	invJ = eye(24,24)*2*(1/cellSize_);
	
	[s, t, p, w] = FEA_GaussianIntegral();
	dShape = FEA_DeShapeFunction(s,t,p);
	eleB = FEA_ElementStrainMatrix(dShape, invJ);
	
	%%Element Elasticity Matrix (D-mat)
	S = FEA_HookeLaw();
	eleD = FEA_ElementElasticityMatrix(S);
	
	%% Integration Weights
	wgt = w.*detJ; wgt = repmat(wgt, 1, 6); wgt = reshape(wgt', 1, numel(wgt));
	Ke = eleB'*(eleD.*wgt)*eleB;
end




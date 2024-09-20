function [compliance, volumeFraction] = FEA_ComputeComplianceVoxel(varargin)
	global meshHierarchy_;
	global F_;
	global U_;
	global strainEnergyPerElement_;
	if isempty(F_), FEA_ApplyBoundaryCondition(); end
	if isempty(meshHierarchy_(1).Ke), FEA_SetupVoxelBased(); end
	% U_ = zeros(meshHierarchy_(1).numDOFs,1);
	tStart = tic;
    if 0==nargin
		densityField = ones(meshHierarchy_(1).numElements,1);		
    else
        densityField = varargin{1};       
    end
	meshHierarchy_(1).eleModulus = TopOpti_MaterialInterpolationSIMP(densityField(:));
	Solving_AssembleFEAstencil();
	Solving_CG_GMGS('printP_ON'); 		
	ceList = TopOpti_ComputeUnitCompliance();
	
	strainEnergyPerElement_ = meshHierarchy_(1).eleModulus(:) .* ceList(:);
	% compliance = meshHierarchy_(1).eleModulus(:)' *ceList(:);
	compliance = sum(strainEnergyPerElement_);
	volumeFraction = sum(densityField)/numel(densityField);
	disp(['Solving Linear System Costs: ' sprintf('%10.1f',toc(tStart)) 's']);
	disp(['Compliance in total (weighted): ' sprintf('%10.5e ', compliance)]);
	disp(['Deposition ratio: ' sprintf('%10.4e ', volumeFraction)]);
end
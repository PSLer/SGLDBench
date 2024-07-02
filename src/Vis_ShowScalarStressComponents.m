function Vis_ShowScalarStressComponents(axHandle, sType, varargin)
	global domainType_;
	global cartesianStressField_;
	global vonMisesStressField_;
	
	switch sType
		case 'Sigma_xx'
			srcField = cartesianStressField_(:,1);
		case 'Sigma_yy'
			srcField = cartesianStressField_(:,2);
		case 'Sigma_zz'
			srcField = cartesianStressField_(:,3);
		case 'Tau_yz'
			srcField = cartesianStressField_(:,4);
		case 'Tau_zx'
			srcField = cartesianStressField_(:,5);
		case 'Tau_xy'
			srcField = cartesianStressField_(:,6);
		case 'Von Mises'
			srcField = vonMisesStressField_;
		otherwise
			warning('Undefined Stress Type!'); return;				
	end	
	if 2==nargin
		Vis_ShowScalarFieldOnVoxelSurface(axHandle, srcField);
	else
		Vis_ShowScalarFieldOnVoxelSurface(axHandle, srcField, varargin{1});
	end	
end
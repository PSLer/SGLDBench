function Vis_ShowDeformation(axHandle, dType, varargin)
	global U_;
	if isempty(U_); return; end
	
	switch dType
		case 'Ux'
			srcField = U_(1:3:end,1);
		case 'Uy'
			srcField = U_(2:3:end,1);
		case 'Uz'
			srcField = U_(3:3:end,1);
		case 'Total'
			srcField = reshape(U_, 3, numel(U_)/3)';
			srcField = vecnorm(srcField,2,2);
		otherwise
			warning('Undefined Deformation Direction!'); return;				
	end	
	if 2==nargin
		Vis_ShowScalarFieldOnVoxelSurface(axHandle, srcField);
	else
		Vis_ShowScalarFieldOnVoxelSurface(axHandle, srcField, varargin{1});
	end	
end
function Vis_ShowScalarFieldOnVoxelSurface(axHandle, scalarField, varargin)
	global meshHierarchy_;
	global U_;
	global boundingBox_;
	
	if 2==nargin
		minFeaterSize = min(boundingBox_(2,:)-boundingBox_(1,:)); selfFac = 10;
		scalingFactor = minFeaterSize/selfFac/max(abs(U_));		
	elseif 3==nargin
		scalingFactor = varargin{1};
	else
		error('Wrong Input!');
	end
	
	deformation = U_ * scalingFactor; %%Exaggerated display
	deformation = reshape(deformation, 3, meshHierarchy_(1).numNodes)';
	deformedCartesianGrid = meshHierarchy_(1).boundaryNodeCoords + deformation(meshHierarchy_(1).nodesOnBoundary,:);
	faces2Draw = meshHierarchy_(1).boundaryEleFaces;
	vertices2Draw = deformedCartesianGrid;
	color2Draw = scalarField(meshHierarchy_(1).nodesOnBoundary,1);
	hd = patch(axHandle, 'Faces', faces2Draw, 'Vertices', vertices2Draw,'FaceVertexCData', color2Draw);

	colormap(axHandle, 'jet');
	set(hd, 'FaceColor', 'interp', 'FaceAlpha', 1, 'EdgeColor', 'none');
	h = colorbar; t=get(h,'Limits'); 
	set(h,'Ticks',linspace(t(1),t(2),7),'AxisLocation','out');	
	L=cellfun(@(x)sprintf('%.2e',x),num2cell(linspace(t(1),t(2),7)),'Un',0); 
	set(h,'xticklabel',L);		
	axis(axHandle, 'off'); axis(axHandle, 'equal'); axis(axHandle, 'tight');
	set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 20);		
end
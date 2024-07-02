function Vis_UserLighting(axHandle)
	% %%Lighting, Reflection
	lighting(axHandle, 'gouraud');
	material(axHandle, 'dull');
	camlight(axHandle, 'headlight', 'infinite');	
end
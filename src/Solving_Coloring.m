function colors = Solving_Coloring(eleMapForward, resX, resY, resZ)
	numColors = 8;
	colors = cell(numColors,1);
	allValidVoxles = reshape(eleMapForward, resY, resX, resZ);
	samp = allValidVoxles(1:2:end,1:2:end,1:2:end); samp = samp(:); colors{1} = samp(samp>0);
	samp = allValidVoxles(2:2:end,1:2:end,1:2:end); samp = samp(:); colors{2} = samp(samp>0);
	samp = allValidVoxles(1:2:end,2:2:end,1:2:end); samp = samp(:); colors{3} = samp(samp>0);
	samp = allValidVoxles(2:2:end,2:2:end,1:2:end); samp = samp(:); colors{4} = samp(samp>0);
	samp = allValidVoxles(1:2:end,1:2:end,2:2:end); samp = samp(:); colors{5} = samp(samp>0);
	samp = allValidVoxles(2:2:end,1:2:end,2:2:end); samp = samp(:); colors{6} = samp(samp>0);
	samp = allValidVoxles(1:2:end,2:2:end,2:2:end); samp = samp(:); colors{7} = samp(samp>0);
	samp = allValidVoxles(2:2:end,2:2:end,2:2:end); samp = samp(:); colors{8} = samp(samp>0);	
end
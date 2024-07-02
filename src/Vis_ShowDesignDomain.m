function Vis_ShowDesignDomain(axHandle)
	global meshHierarchy_;
	global passiveElements_;
	
	FV.faces = meshHierarchy_(1).boundaryEleFaces;
	FV.vertices = meshHierarchy_(1).boundaryNodeCoords;
	hdSilhouette = patch(axHandle, FV);
	refVertices = [
		-0.5	-0.5	-0.5
		0.5		-0.5	-0.5
		0.5		0.5		-0.5
		-0.5	0.5		-0.5
		-0.5	-0.5	0.5
		0.5		-0.5	0.5
		0.5		0.5		0.5
		-0.5	0.5		0.5		
	];	
	
	if ~isempty(passiveElements_)
		drawPassiveEles.vertices = zeros(numel(passiveElements_)*8,3);
		drawPassiveEles.faces = zeros(numel(passiveElements_), 8);
		for ii=1:numel(passiveElements_)
			for jj=1:8
				drawPassiveEles.vertices((ii-1)*8+jj,:) = refVertices(jj,:) + meshHierarchy_(1).eleCentroidList(passiveElements_(ii),:);
			end
			drawPassiveEles.faces(ii,:) = (1:8) + 8*(ii-1);
		end
		drawPassiveEles.faces = drawPassiveEles.faces(:,[4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
		drawPassiveEles.faces = reshape(drawPassiveEles.faces(:), 4, 6*numel(passiveElements_))';	
		hold('on');
		hdPassive = patch(axHandle, drawPassiveEles); 
	else
		hdPassive = [];
	end
	set(hdSilhouette, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'None');
	set(hdPassive, 'FaceColor', [252 141 98]/255, 'FaceAlpha', 1.0, 'EdgeColor', 'None');

	% lighting('gouraud');
	% material('dull'); %%shiny; dull; metal	
	% camlight('headlight','infinite');
	
	axis(axHandle, 'equal'); axis(axHandle, 'tight'); axis(axHandle, 'on');
	% set(axHandle, 'FontName', 'Times New Roman', 'FontSize', 20);		
end
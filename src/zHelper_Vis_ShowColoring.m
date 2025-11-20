function zHelper_Vis_ShowColoring(iLevel, varargin)
	global meshHierarchy_;
	
	nx = meshHierarchy_(1).resX;
	ny = meshHierarchy_(1).resY;
	nz = meshHierarchy_(1).resZ;
	numElements = meshHierarchy_(iLevel).numElements;
	numNodes = meshHierarchy_(iLevel).numNodes;
	eNodMat = meshHierarchy_(iLevel).eNodMat;
	
	numAllNodes = (1+nx)*(1+ny)*(1+nz);
	xSeed = 0:nx; xSeed = single(xSeed);
	ySeed = ny:-1:0; ySeed = single(ySeed);
	zSeed = 0:nz; zSeed = single(zSeed);
	nodeCoords = zeros(numNodes,3, 'single');
	nodMapBack = meshHierarchy_(iLevel).nodMapBack;
	tmp = repmat(reshape(repmat(xSeed,ny+1,1), (nx+1)*(ny+1), 1), (nz+1), 1); 
	nodeCoords(:,1) = tmp(nodMapBack,1);
	tmp = repmat(repmat(ySeed,1,nx+1)', (nz+1), 1);
	nodeCoords(:,2) = tmp(nodMapBack,1);
	tmp = reshape(repmat(zSeed,(nx+1)*(ny+1),1), (nx+1)*(ny+1)*(nz+1), 1);
	nodeCoords(:,3) = tmp(nodMapBack,1);
	
	numNod2ElesVec = zeros(numNodes,1,'int32');
	for jj=1:8
		iNodes = eNodMat(:,jj);
		numNod2ElesVec(iNodes) = numNod2ElesVec(iNodes) + 1;
	end	
	nodeState = zeros(numNodes,1,'int32');
	nodeState(numNod2ElesVec<8) = 1;	
	tmp = zeros(numElements,1,'int32');
	for ii=1:8
		tmp = tmp + nodeState(eNodMat(:,ii));
	end
	elementsOnBoundary = int32(find(tmp>0));
	
	if 1==nargin
		numColors = numel(meshHierarchy_(iLevel).colors); % ==8
		colors_faces = cell(numColors,1);
		for ii=1:numColors
			iColorBoundaryElements = intersect(elementsOnBoundary, meshHierarchy_(iLevel).colors{ii});
			iPatchIndices = eNodMat(iColorBoundaryElements, [4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
			iPatchIndices = reshape(iPatchIndices(:), 4, 6*numel(iColorBoundaryElements));
			tmp2 = nodeState(iPatchIndices); tmp2 = sum(tmp2,1);
			iBoundaryEleFaces = iPatchIndices(:,4==tmp2);
			colors_faces{ii} = iBoundaryEleFaces';
		end
		
		figure;
		iColorPatch.vertices = nodeCoords;
		
		iColorPatch.faces = colors_faces{1};
		hd_1 = patch(iColorPatch); hold on;
		
		iColorPatch.faces = colors_faces{2};
		hd_2 = patch(iColorPatch); hold on;	
	
		iColorPatch.faces = colors_faces{3};
		hd_3 = patch(iColorPatch); hold on;	
		
		iColorPatch.faces = colors_faces{4};
		hd_4 = patch(iColorPatch); hold on;	
		
		iColorPatch.faces = colors_faces{5};
		hd_5 = patch(iColorPatch); hold on;	
		
		iColorPatch.faces = colors_faces{6};
		hd_6 = patch(iColorPatch); hold on;	
		
		iColorPatch.faces = colors_faces{7};
		hd_7 = patch(iColorPatch); hold on;	
		
		iColorPatch.faces = colors_faces{8};
		hd_8 = patch(iColorPatch); hold on;	
		
		set(hd_1, 'FaceColor', [227 26 28]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		set(hd_2, 'FaceColor', [31 120 180]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		set(hd_3, 'FaceColor', [251 154 153]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		set(hd_4, 'FaceColor', [166 206 227]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		set(hd_5, 'FaceColor', [255 127 0]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		set(hd_6, 'FaceColor', [51 160 44]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		set(hd_7, 'FaceColor', [106 61 154]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		set(hd_8, 'FaceColor', [178 223 138]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		% set([hd_1; hd_2; hd_3; hd_4; hd_5; hd_6; hd_7; hd_8;], 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
	else
		numColors = numel(meshHierarchy_(iLevel).colors); % ==8
		colors_faces = cell(numColors,1);
		allElements = int32(1:numElements)';
		iColor = varargin{1};
		iColorBoundaryElements = intersect(allElements, meshHierarchy_(iLevel).colors{iColor});
		iPatchIndices = eNodMat(iColorBoundaryElements, [4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
		iPatchIndices = reshape(iPatchIndices(:), 4, 6*numel(iColorBoundaryElements));
		colors_faces = iPatchIndices';	
		
		figure;
		iColorPatch.vertices = nodeCoords;
		iColorPatch.faces = colors_faces;
		hd = patch(iColorPatch); hold on;
		%%1: [227 26 28]/255
		%%2: [31 120 180]/255
		switch iColor
			case 1
				set(hd, 'FaceColor', [227 26 28]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
			case 2
				set(hd, 'FaceColor', [31 120 180]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
			case 3
				set(hd, 'FaceColor', [251 154 153]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
			case 4
				set(hd, 'FaceColor', [166 206 227]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
			case 5
				set(hd, 'FaceColor', [255 127 0]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
			case 6
				set(hd, 'FaceColor', [51 160 44]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
			case 7
				set(hd, 'FaceColor', [106 61 154]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
			case 8
				set(hd, 'FaceColor', [178 223 138]/255, 'FaceAlpha', 1, 'EdgeColor', 'k');
		end
		
	end
	
	view(35,10);
	camlight('headlight');
	axis('equal'); axis('tight'); axis('off');			
end
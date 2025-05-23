function Vis_ShowPSLsSmoothed(axHandle)
	global meshHierarchy_;
	global boundingBoxSmoothed_;
	global silhouetteStructSmoothed_;
	global majorPSLpoolSmoothed_;
	global mediumPSLpoolSmoothed_;
	global minorPSLpoolSmoothed_;
	% global silhouetteStruct_;
	
	numTarMajorPSLs = numel(majorPSLpoolSmoothed_);
	numTarMediumPSLs = numel(mediumPSLpoolSmoothed_);
	numTarMinorPSLs = numel(minorPSLpoolSmoothed_);
	color4MajorPSLs = struct('arr', []); color4MajorPSLs = repmat(color4MajorPSLs, numTarMajorPSLs, 1);
	color4MediumPSLs = struct('arr', []); color4MediumPSLs = repmat(color4MediumPSLs, numTarMediumPSLs, 1);
	color4MinorPSLs = struct('arr', []); color4MinorPSLs = repmat(color4MinorPSLs, numTarMinorPSLs, 1);	
	
	for ii=1:numTarMajorPSLs
		color4MajorPSLs(ii).arr = ones(1, majorPSLpoolSmoothed_(ii).length);
	end
	for ii=1:numTarMediumPSLs
		color4MediumPSLs(ii).arr = ones(1, mediumPSLpoolSmoothed_(ii).length);
	end			
	for ii=1:numTarMinorPSLs
		color4MinorPSLs(ii).arr = ones(1, minorPSLpoolSmoothed_(ii).length);
	end	

	thicknessScaling = 100;
	lineWidthTube = min(boundingBoxSmoothed_(2,:)-boundingBoxSmoothed_(1,:))/thicknessScaling;
	[gridXmajor, gridYmajor, gridZmajor, gridCmajor, ~] = ExpandPSLs2Tubes(majorPSLpoolSmoothed_, color4MajorPSLs, lineWidthTube);
	[gridXmedium, gridYmedium, gridZmedium, gridCmedium, ~] = ExpandPSLs2Tubes(mediumPSLpoolSmoothed_, color4MediumPSLs, lineWidthTube);
	[gridXminor, gridYminor, gridZminor, gridCminor, ~] = ExpandPSLs2Tubes(minorPSLpoolSmoothed_, color4MinorPSLs, lineWidthTube);

	%%Show silhouette
	handleMajorPSL = []; handleMediumPSL = []; handleMinorPSL = [];
	hSilo = patch(gca, silhouetteStructSmoothed_); hold(gca, 'on');	
	if ~isempty(gridXmajor)
		hold(gca, 'on'); 
		handleMajorPSL = surf(axHandle, gridXmajor, gridYmajor, gridZmajor, gridCmajor);
	end
	if ~isempty(gridXmedium)
		hold(gca, 'on');
		handleMediumPSL = surf(axHandle, gridXmedium, gridYmedium, gridZmedium, gridCmedium);
	end
	if ~isempty(gridXminor)
		hold(gca, 'on');
		handleMinorPSL = surf(axHandle, gridXminor, gridYminor, gridZminor, gridCminor);
	end
	set(hSilo, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'None');
	set(handleMajorPSL, 'FaceColor', [191 129 45]/255, 'EdgeColor', 'None');
	set(handleMediumPSL, 'FaceColor', [148 0 211]/255, 'EdgeColor', 'None');
	set(handleMinorPSL, 'FaceColor', [53 151 143]/255, 'EdgeColor', 'None');	
	set(hSilo, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'None');
	axis(axHandle, 'equal'); axis(axHandle, 'tight'); axis(axHandle, 'off');
end

function [gridX, gridY, gridZ, gridC, gridIndices] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	%%Syntax
	%% [gridX, gridY, gridZ, gridC] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	gridX = [];
	gridY = [];
	gridZ = [];
	gridC = [];
	gridIndices = [];
	if isempty(PSLs), return; end
	n = 8; 
	numLines = length(PSLs);
	gridXYZ = zeros(3,n+1,1);
	gridC = zeros(n+1,1);
	gridIndices = struct('arr', []);
	gridIndices = repmat(gridIndices, numLines, 1);	
	for ii=1:numLines		
		curve = PSLs(ii).phyCoordList';
		npoints = size(curve,2);
		%deltavecs: average for internal points. first strecth for endpoitns.		
		dv = curve(:,[2:end,end])-curve(:,[1,1:end-1]);		
		%make nvec not parallel to dv(:,1)
		nvec=zeros(3,1); 
		[~,idx]=min(abs(dv(:,1))); 
		nvec(idx)=1;
		%precalculate cos and sing factors:
		cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
		sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
		%Main loop: propagate the normal (nvec) along the tube
		xyz = zeros(3,n+1,npoints+2);
		for k=1:npoints
			convec=cross(nvec,dv(:,k));
			convec=convec./norm(convec);
			nvec=cross(dv(:,k),convec);
			nvec=nvec./norm(nvec);
			%update xyz:
			xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1]) + cfact.*repmat(r*nvec,[1,n+1]) + sfact.*repmat(r*convec,[1,n+1]);
        end
		%finally, cap the ends:
		xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
		xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
		gridIndices(ii).arr = size(gridXYZ,3) : size(gridXYZ,3)+size(xyz,3)-1;
		gridXYZ(:,:,end+1:end+npoints+2) = xyz;	
		color = colorSrc(ii).arr;	
		c = [color(1) color color(end)];
		c = repmat(c, n+1, 1);
		gridC(:,end+1:end+npoints+2) = c;
	end		
	gridX = squeeze(gridXYZ(1,:,:)); 
	gridX(:,1) = [];
	gridY = squeeze(gridXYZ(2,:,:)); 
	gridY(:,1) = [];
	gridZ = squeeze(gridXYZ(3,:,:)); 
	gridZ(:,1) = [];
	gridC(:,1) = [];
end
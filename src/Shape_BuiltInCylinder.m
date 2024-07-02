function Shape_BuiltInCylinder(radius, ratioInnerRadius, height)
	global surfaceTriMesh_;
	[vertices, faces] = CreateCylinder(radius, ratioInnerRadius, height);
	surfaceTriMesh_.nodeCoords = vertices;
	surfaceTriMesh_.eNodMat = faces;
	surfaceTriMesh_.numNodes = size(surfaceTriMesh_.nodeCoords,1); 
	surfaceTriMesh_.numElements = size(surfaceTriMesh_.eNodMat,1);
end


function [vertices, faces] = CreateCylinder(radius, ratioInnerRadius, height)
	nnp = 64;
	
	vertices_plate_outer_bottom = zeros(nnp, 3);
	for ii=1:nnp
		iRad = ii/nnp * 2*pi;
		vertices_plate_outer_bottom(ii,1:2) = radius*[cos(iRad) sin(iRad)];
	end
	vertices_plate_outer_top = vertices_plate_outer_bottom;
	vertices_plate_outer_top(:,3) = height;
	
	if 0==ratioInnerRadius
		ctr_bottom = [0 0 0];
		ctr_top = [0 0 height];
		vertices = [vertices_plate_outer_bottom; vertices_plate_outer_top; ctr_bottom; ctr_top];
		faces = zeros(4*nnp,3);
		faces(1:nnp,1) = 2*nnp + 1;
		faces(1:nnp,2) = [(2:nnp) 1]';
		faces(1:nnp,3) = (1:nnp)';

		faces(nnp+(1:nnp),1) = 2*nnp + 2;
		faces(nnp+(1:nnp),2) = nnp+(1:nnp)';	
		faces(nnp+(1:nnp),3) = (nnp+[2:nnp, 1])';
		
		tmp = [(1:nnp)', [(2:nnp)'; 1], nnp+[(2:nnp)'; 1], nnp+(1:nnp)'];
		tmp = tmp(:,[1 2 3 3 4 1])'; tmp = tmp(:)'; tmp = reshape(tmp,3,2*nnp)';
		faces((2*nnp+1):4*nnp,:) = tmp;
	else
		vertices_plate_inner_bottom = vertices_plate_outer_bottom * ratioInnerRadius;
		vertices_plate_inner_top = vertices_plate_inner_bottom;
		vertices_plate_inner_top(:,3) = height;
		vertices = [vertices_plate_outer_bottom; vertices_plate_outer_top; vertices_plate_inner_bottom; vertices_plate_inner_top];

		tmp = [(1:nnp)', [(2:nnp)'; 1], nnp+[(2:nnp)'; 1], nnp+(1:nnp)'];
		tmp = tmp(:,[1 2 3 3 4 1])'; tmp = tmp(:)'; tmp = reshape(tmp,3,2*nnp)';
		faces = tmp;
		tmp = [(1:nnp)', nnp+(1:nnp)', nnp+[(2:nnp)'; 1], [(2:nnp)'; 1]] + 2*nnp;
		tmp = tmp(:,[1 2 3 3 4 1])'; tmp = tmp(:)'; tmp = reshape(tmp,3,2*nnp)';
		faces = [faces; tmp];
		tmp = [(1:nnp)', 2*nnp+(1:nnp)', 2*nnp+[2:nnp, 1]', [2:nnp, 1]'];
		tmp = tmp(:,[1 2 3 3 4 1])'; tmp = tmp(:)'; tmp = reshape(tmp,3,2*nnp)';
		faces = [faces; tmp];
		tmp = [nnp+(1:nnp)', nnp+[2:nnp, 1]', 3*nnp+[2:nnp, 1]', 3*nnp+(1:nnp)'];
		tmp = tmp(:,[1 2 3 3 4 1])'; tmp = tmp(:)'; tmp = reshape(tmp,3,2*nnp)';
		faces = [faces; tmp];
	end
end
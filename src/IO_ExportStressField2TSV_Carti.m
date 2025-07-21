function IO_ExportStressField2TSV_Carti(fileName, cartesianStresses, xPhys, threshold)
	global meshHierarchy_;
	global boundingBox_;
	global cellSize_; 
	solidEles = find(xPhys>threshold);
	numElementsNew = numel(solidEles);
	eNodMatNew = meshHierarchy_(1).eNodMat(solidEles,:);
	newNodes = unique(eNodMatNew);
	numNodesNew = numel(newNodes);
	nodMap = zeros(meshHierarchy_(1).numNodes,1);
	nodMap(newNodes,1) = (1:numNodesNew)';
	eNodMatNew = nodMap(eNodMatNew);
	cartesianStressesNew = cartesianStresses(newNodes,:);
	
	%fid = fopen(strcat(outPath_, fileName, '.TSVcarti'), 'w');
	fid = fopen(fileName, 'w');
	%%Write in Economic Cartesian Mesh Style
	%%2.1 File Header
	fprintf(fid, '%s %s', 'DATASET CARTESIAN_GRID'); fprintf(fid, '\n');
    fprintf(fid, '%s %s %s %s', 'Stress Data Type: NODE'); fprintf(fid, '\n');
	fprintf(fid, '%s', 'Resolution:');
	disp('Writing Stress Field to .TSVcarti ......');
	%%2.2 Mesh Description
	fprintf(fid, ' %d %d %d\n', [meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]);
	fprintf(fid, '%s', 'LowerBound:');
	fprintf(fid, ' %.6f %.6f %.6f\n', boundingBox_(1,:));
	fprintf(fid, '%s', 'UpperBound:');
	fprintf(fid, ' %.6f %.6f %.6f\n', boundingBox_(2,:)*cellSize_);	
	fprintf(fid, '%s', 'ELEMENTS');
	fprintf(fid, ' %d', numElementsNew);
	fprintf(fid, ' %s \n', 'int');	
	fprintf(fid, '%d\n', meshHierarchy_(1).eleMapBack(solidEles,1)');
	%%2.2 Cartesian Stress 
	fprintf(fid, '%s %s %s %s ', 'Number of Stress Fields:');
	fprintf(fid, '%d\n', 1);
	fprintf(fid, '%s %s ', 'Node Forces:'); 
	fprintf(fid, '%d\n', 0);
	fprintf(fid, '%s %s ', 'Fixed Nodes:'); fprintf(fid, '%d\n', 0);
	fprintf(fid, '%s %s ', 'Cartesian Stress:'); fprintf(fid, '%d\n', numNodesNew);
	fprintf(fid, '%.6e %.6e %.6e %.6e %.6e %.6e\n', cartesianStressesNew');	
	fclose(fid);
	disp('Done with Writing Stress Field to .TSVcarti!');
end
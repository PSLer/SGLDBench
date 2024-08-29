function IO_WriteGraph_OBJ(fileName, nodeCoords, edges)

	fid = fopen(fileName, 'w');
	for ii=1:size(nodeCoords,1)
		fprintf(fid, '%s ', 'v');
		fprintf(fid, '%.6f %.6f %.6f\n', nodeCoords(ii,:));
	end
	for ii=1:size(edges,1)
		fprintf(fid, '%s ', 'l');
		fprintf(fid, '%d %d\n', edges(ii,:));
	end		
	fclose(fid);	
end
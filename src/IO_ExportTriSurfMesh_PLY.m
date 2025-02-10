function IO_ExportTriSurfMesh_PLY(nodeCoords, eNodMat, fileName)
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s\n', 'ply');
	fprintf(fid, '%s %s', 'format ascii'); fprintf(fid, '%.1f\n', 1.0);
	fprintf(fid, '%s %s', 'element vertex'); fprintf(fid, ' %d\n', size(nodeCoords,1));
	fprintf(fid, '%s %s %s', 'property float x'); fprintf(fid, '\n');
	fprintf(fid, '%s %s %s', 'property float y'); fprintf(fid, '\n'); 
	fprintf(fid, '%s %s %s', 'property float z'); fprintf(fid, '\n');
	fprintf(fid, '%s %s', 'element face'); fprintf(fid, ' %d\n', size(eNodMat,1));
	fprintf(fid, '%s %s %s %s %s', 'property list uchar uint vertex_indices'); fprintf(fid, '\n');
	fprintf(fid, '%s\n', 'end_header');
	fprintf(fid, '%.6f %.6f %.6f\n', nodeCoords');
	fprintf(fid, '%d %d %d %d\n', [repmat(3,size(eNodMat,1),1) eNodMat-1]');					
	fclose(fid);	
end
function IO_ExportStressField2TSV()
	global outPath_;
	global stressFieldOnTetMesh_; 
	global gateWayTetMesh_;	
	
	fileName = strcat(outPath_, 'StressField_Tet_v2.stress');
	fid = fopen(fileName, 'w');
	fprintf(fid, '%s ', 'Version');
	fprintf(fid, '%.1f\n', 2.0);
	
	fprintf(fid, '%s %s ', 'Solid Tet');
	fprintf(fid, '%d\n', 1);
	
	fprintf(fid, '%s ', 'Vertices:');
	fprintf(fid, '%d\n', gateWayTetMesh_.numNodes);		
	fprintf(fid, '%.6e %.6e %.6e\n', gateWayTetMesh_.nodeCoords');

	fprintf(fid, '%s ', 'Elements:');
	fprintf(fid, '%d \n', gateWayTetMesh_.numElements);
	fprintf(fid, '%d %d %d %d\n', gateWayTetMesh_.eNodMat');

	fprintf(fid, '%s %s ', 'Node Forces:'); 
	fprintf(fid, '%d\n', 0);
	fprintf(fid, '%s %s ', 'Fixed Nodes:'); fprintf(fid, '%d\n', 0);

	fprintf(fid, '%s %s', 'Cartesian Stress:'); 
	fprintf(fid, '%d\n', gateWayTetMesh_.numNodes);
	fprintf(fid, '%.6e %.6e %.6e %.6e %.6e %.6e\n', stressFieldOnTetMesh_');
	fclose(fid);
end
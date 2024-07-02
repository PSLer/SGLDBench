function IO_ExportTopOptiResults()
	global outPath_;
	global meshHierarchy_;
	global loadingCond_; global fixingCond_; 
	global cHist_; global volHist_; global sharpHist_; global densityLayout_;

	%%1. Compliance
	fid = fopen(strcat(outPath_, 'complianceHist.dat'), 'w');
	fprintf(fid, '%16.6e\n', cHist_);
	fclose(fid);
	
	%%2. Volume Fraction
	fid = fopen(strcat(outPath_, 'volumeFractionHist.dat'), 'w');
	fprintf(fid, '%16.6e\n', volHist_);
	fclose(fid);

	%%3. sharpness
	fid = fopen(strcat(outPath_, 'sharpnessHist.dat'), 'w');
	fprintf(fid, '%16.6e\n', sharpHist_);
	fclose(fid);	
		
	%%4. optimized result in SIMP for standard FEA
	fid = fopen(strcat(outPath_, 'optimizedModel.topopti'), 'w');
	fprintf(fid, '%s %s %s', '# topology optimization'); fprintf(fid, '\n');
	fprintf(fid, '%s %s', 'domain type:');
	fprintf(fid, ' %s', '3D'); fprintf(fid, '\n');
	fprintf(fid, '%s', 'resolution:');
	fprintf(fid, ' %d %d %d\n', [meshHierarchy_(1).resX meshHierarchy_(1).resY meshHierarchy_(1).resZ]);
	fprintf(fid, '%s %s %s', 'optimized material density:');
	fprintf(fid, ' %d\n', meshHierarchy_(1).numElements);
	fprintf(fid, '%d %.6f\n', [double(meshHierarchy_(1).eleMapBack) densityLayout_]');
	fprintf(fid, '%s %s', 'fixed position:');
	fprintf(fid, ' %d\n', size(fixingCond_,1));
	fprintf(fid, '%d %d %d %d\n', [meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary(fixingCond_(:,1))) fixingCond_(:,2:end)]');
	fprintf(fid, '%s %s', 'loading condition:');
	fprintf(fid, ' %d\n', size(loadingCond_,1));
	fprintf(fid, '%d %.2e %.2e %.2e\n', [double(meshHierarchy_(1).nodMapBack(meshHierarchy_(1).nodesOnBoundary(loadingCond_(:,1)))) loadingCond_(:,2:end)]');
	fprintf(fid, '%s %s %s', 'additional boundary conditions:');
	fprintf(fid, ' %d\n', 0);
	fclose(fid);
end
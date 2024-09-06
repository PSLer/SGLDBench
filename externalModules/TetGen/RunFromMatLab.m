
maxCellVol = 10.123;
callTetGen = strcat('"./tetgen.exe" -gq1.414a', num2str(maxCellVol), ' gatewayMesh.ply');
system(callTetGen);
function IO_ExtractReconstructExportDesignIsosurface(oFileName)
	global meshHierarchy_;
	global densityLayout_;
	global inputSolidMesh_;
	global surfaceTriMesh_;
	global FVreconstructed_;
	global outPath_;
	
	%%Extract Isosurface
	disp('Extracting Isosurface of Design ...');
	solidEles = find(densityLayout_>0.1);
	exteriorEles = meshHierarchy_(1).eleMapBack(solidEles);
	valForIIsosurface = zeros(meshHierarchy_(1).resX*meshHierarchy_(1).resY*meshHierarchy_(1).resZ,1);
	valForIIsosurface(exteriorEles) = 1;
	valForIIsosurface = reshape(valForIIsosurface, meshHierarchy_(1).resY, meshHierarchy_(1).resX, meshHierarchy_(1).resZ);
	valForIIsosurface = flip(valForIIsosurface,1);
	valForIIsosurface = smooth3(valForIIsosurface,'box',1);
	iSurface = isosurface(valForIIsosurface,0);
	iCap = isocaps(valForIIsosurface,0);	
	iCap.faces = size(iSurface.vertices,1) + iCap.faces;
	FV.vertices = [iSurface.vertices; iCap.vertices];
	FV.faces = [iSurface.faces; iCap.faces];				
	FV = CompactPatchVertices(FV);		
	
	%%Tet-meshing TopOpti Design
	dataName = strcat(outPath_, 'TopOptiDesignIsosurface_RAW.ply');
	IO_ExportTriSurfMesh_PLY(FV.vertices, FV.faces, dataName);
	if exist('../externalModules/TetGen/', 'dir')
		callTetWild_Executable = strcat('"../externalModules/fTetWild/FloatTetwild_bin.exe" -i', ...
			char(strcat(" ", dataName, " ", '--no-binary', " ", '-l', " ", '0.05'))); 
	else
		callTetWild_Executable = strcat('"./externalModules/fTetWild/FloatTetwild_bin.exe" -i', ...
			char(strcat(" ", dataName, " ", '--no-binary', " ", '-l', " ", '0.05'))); 	
	end
	system(callTetWild_Executable); pause(1);
	
	%%Extracting Reconstructed Surface
	odataName = strcat(outPath_, 'TopOptiDesignIsosurface_RAW.ply_.msh');
	backupSurfaceMesh = surfaceTriMesh_;
	backSolidMesh = inputSolidMesh_;
	IO_ImportSolidMesh(odataName);
	inputSolidMesh_ = backSolidMesh;
	FVreconstructed_ = surfaceTriMesh_;
	surfaceTriMesh_ = backupSurfaceMesh;
	IO_ExportTriSurfMesh_PLY(FVreconstructed_.nodeCoords, FVreconstructed_.eNodMat, oFileName);
	disp('...Done with Isosurface Re-construction!');
end

function oPatchs = CompactPatchVertices(iPatchs)
	oPatchs = iPatchs;
	numOriginalVertices = size(iPatchs.vertices,1);
	numOriginalFaces = size(iPatchs.faces,1);
	validVertices = unique(iPatchs.faces);
	numValidVertices = size(validVertices,1);
	if numOriginalFaces==numOriginalVertices, return; end
	mapVerticesValid2Original = zeros(numOriginalVertices,1);
	mapVerticesValid2Original(validVertices) = (1:numValidVertices)';
	oPatchs.vertices = oPatchs.vertices(validVertices,:);
	oPatchs.faces = mapVerticesValid2Original(iPatchs.faces);
end
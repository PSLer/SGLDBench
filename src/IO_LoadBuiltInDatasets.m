function IO_LoadBuiltInDatasets(MdlSelect)
	global loadingCond_;
	global fixingCond_;
	%% !!!
	%%The boundary condition files "x.bc" are created by directly copying the values of the global variables "loadingCond_" and "fixingCond_"
	%% !!!
	switch MdlSelect
		case 'Bone'
			IO_ImportSurfaceMesh('../data/Bone.ply');
			if 1 %% Target Resolution: 512
				FEA_CreateVoxelizedModel(512);
				FEA_VoxelBasedDiscretization();
				loadingCond_ = load('../data/Bone_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
				fixingCond_ = load('../data/Bone_R512_fixa.bc');
			else %% Target Resolution: 1200
				FEA_CreateVoxelizedModel(1200);
				FEA_VoxelBasedDiscretization();
				loadingCond_ = load('../data/Bone_R1200_loads.bc'); %%Load prescribed boundary conditions for TESTING
				fixingCond_ = load('../data/Bone_R1200_fixa.bc');		
			end
		case 'Part'
			if 1
				IO_ImportSurfaceMesh('../data/Part.ply');
				FEA_CreateVoxelizedModel(512);
				FEA_VoxelBasedDiscretization();
				loadingCond_ = load('../data/Part_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
				fixingCond_ = load('../data/Part_R512_fixa.bc');
			else
				IO_ImportTopVoxels('../data/Part_R256.TopVoxel'); %%Create from wrapped voxel file
			end
		case 'Part2'
			IO_ImportSurfaceMesh('../data/Part2.ply');
			FEA_CreateVoxelizedModel(512);
			FEA_VoxelBasedDiscretization();
			loadingCond_ = load('../data/Part2_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
			fixingCond_ = load('../data/Part2_R512_fixa.bc');
		case 'Part3'
			IO_ImportSurfaceMesh('../data/Part3.ply');
			FEA_CreateVoxelizedModel(512);
			FEA_VoxelBasedDiscretization();
			loadingCond_ = load('../data/Part3_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
			fixingCond_ = load('../data/Part3_R512_fixa.bc');		
		case 'Bracket_GE'
			IO_ImportSurfaceMesh('../data/bracket_GE.ply');
			FEA_CreateVoxelizedModel(512);
			FEA_VoxelBasedDiscretization();
			loadingCond_ = load('../data/bracket_GE_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
			fixingCond_ = load('../data/bracket_GE_R512_fixa.bc');			
		case 'Molar'
			IO_ImportSurfaceMesh('../data/Molar.ply');
			FEA_CreateVoxelizedModel(512);
			FEA_VoxelBasedDiscretization();
			loadingCond_ = load('../data/Molar_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
			fixingCond_ = load('../data/Molar_R512_fixa.bc');	
		case 'Fertility'
			IO_ImportSurfaceMesh('../data/Fertility.ply');
			FEA_CreateVoxelizedModel(512);
			FEA_VoxelBasedDiscretization();
			loadingCond_ = load('../data/Fertility_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
			fixingCond_ = load('../data/Fertility_R512_fixa.bc');		
		case 'Hanger'
			IO_ImportSurfaceMesh('../data/Hanger.ply');
			FEA_CreateVoxelizedModel(512);
			FEA_VoxelBasedDiscretization();
			loadingCond_ = load('../data/Hanger_R512_loads.bc'); %%Load prescribed boundary conditions for TESTING
			fixingCond_ = load('../data/Hanger_R512_fixa.bc');
		case 'TopOptiShape'
			IO_ImportSurfaceMesh('../data/TopOptiShape.ply');
			FEA_CreateVoxelizedModel(600);
			FEA_VoxelBasedDiscretization();
			loadingCond_ = load('../data/TopOptiShape_R600_loads.bc'); %%Load prescribed boundary conditions for TESTING
			fixingCond_ = load('../data/TopOptiShape_R600_fixa.bc');		
		otherwise
			error('Please Specify your Dataset!');
	end
end
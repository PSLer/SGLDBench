function SAGS_ComputingStressFieldOnTetMesh()
	global stressFieldOnTetMesh_; 
	global gateWayTetMesh_;
	global dataPrep4SAGS_;
	global modulus_;
    global poissonRatio_;
    global complianceSolid_;
    
	solidMesh4Sim = gateWayTetMesh_;
	%FEA with pure HEX or TET mesh
	[loads, fixations] = Common_TransferBCsFromVoxels2ArbitraryMesh(solidMesh4Sim);
	[deformation, stressFieldOnTetMesh_, compliance] = FEA_SolidTetOrHexMesh(solidMesh4Sim, loads, fixations, modulus_, poissonRatio_);
	if 1
		disp(['Compliance (Voxel): ', sprintf('%f', complianceSolid_)]);
		disp(['Compliance (Hex/Tet): ', sprintf('%f', compliance)]);	
	end	
	
	dataPrep4SAGS_ = struct('nodeCoords', [], 'eNodMat', [], 'cartesianStress', [], 'ps', [], 'frameField', []);
	dataPrep4SAGS_.nodeCoords = gateWayTetMesh_.nodeCoords;
	dataPrep4SAGS_.eNodMat = gateWayTetMesh_.eNodMat;
	dataPrep4SAGS_.cartesianStress = stressFieldOnTetMesh_;
	dataPrep4SAGS_.ps = zeros(size(dataPrep4SAGS_.nodeCoords,1),12);
	for ii=1:size(dataPrep4SAGS_.nodeCoords,1)
		dataPrep4SAGS_.ps(ii,:)= FEA_ComputePrincipalStress(dataPrep4SAGS_.cartesianStress(ii,:));
	end	
	
	%%Write Gateway Stress File for Graded Voronoi Diagram Generation
	IO_ExportStressField2TSV();
end


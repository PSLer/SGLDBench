function [deformation, cartesianStress, compliance] = FEA_SolidTetOrHexMesh(solidMesh4Sim, loads, fixations, E0, v0)
	
	%%Apply for Boundary Conditions
	numDOFs_ = solidMesh4Sim.numNodes * 3;
	freeDOFs_ = (1:numDOFs_)';
	tmp = 3*fixations(:,1);
	fixedDOFs = [tmp-2 tmp-1 tmp]'; fixedDOFs = fixedDOFs(:);
	fixingState = fixations(:,2:end)';
	realFixedDOFs = fixedDOFs(find(1==fixingState(:)));
	freeDOFs_ = setdiff(freeDOFs_, realFixedDOFs);
	Fvec = sparse(numDOFs_,1);
	deformation = zeros(size(freeDOFs_));
	tmp = 3*loads(:,1);
	loadedDOFs_ = [tmp-2 tmp-1 tmp]'; loadedDOFs_ = loadedDOFs_(:);
	loadingVec = loads(:,2:4)'; loadingVec = loadingVec(:);
	Fvec(loadedDOFs_) = loadingVec;
	Fvec = Fvec(freeDOFs_);
	
	switch solidMesh4Sim.meshType
		case 'HEX'
			nEND = 3;
			nEGIP = 8;
            nEleNodes = 8;
			tmp = 3*solidMesh4Sim.eNodMat; eDofMat_ = [tmp-2 tmp-1 tmp];
			eDofMat_ = eDofMat_(:,[1 9 17 2 10 18 3 11 19 4 12 20 5 13 21 6 14 22 7 15 23 8 16 24]);				
		case 'TET'
			nEND = 3;
			nEGIP = 4;
        	nEleNodes = 4;
			tmp = 3*solidMesh4Sim.eNodMat; eDofMat_ = [tmp-2 tmp-1 tmp];
			eDofMat_ = eDofMat_(:,[1 5 9  2 6 10  3 7 11  4 8 12]);			
	end
	
	%%Assemble Stiffness Matrix
	intPots = GaussianIntegral(solidMesh4Sim.meshType);
	gaussIPs = intPots(1:3,:)';
	wgts = intPots(4,:)';
	deShapeFuncs_ = DeShapeFunction(gaussIPs, solidMesh4Sim.meshType);	
	iInvJ = struct('arr', sparse(nEND*nEGIP,nEND*nEGIP));
	iDetJ = zeros(nEGIP,1);
	invJ_ = repmat(iInvJ, solidMesh4Sim.numElements, 1);
	detJ_ = repmat(iDetJ,1,solidMesh4Sim.numElements);
	for ii=1:solidMesh4Sim.numElements
		probeEleNods = solidMesh4Sim.nodeCoords(solidMesh4Sim.eNodMat(ii,:)',:);
		for kk=1:nEGIP
			Jac = deShapeFuncs_(nEND*(kk-1)+1:nEND*kk,:)*probeEleNods;
			iInvJ.arr(nEND*(kk-1)+1:nEND*kk, nEND*(kk-1)+1:nEND*kk) = inv(Jac);
			iDetJ(kk) = det(Jac);	
		end
		invJ_(ii) = iInvJ;
		detJ_(:,ii) = iDetJ;
	end		
	
	[eKi, eKj, eKk] = GetLowerEleStiffMatIndices(solidMesh4Sim.meshType);
	numEntries = length(eKk);
	eleDmat = ElementElasticityMatrix(E0, v0, solidMesh4Sim.meshType);
	sK = zeros(numEntries, solidMesh4Sim.numElements);
	index = 0;
	for ii=1:solidMesh4Sim.numElements
		index = index + 1;						
		eleBmat = ElementStrainMatrix(deShapeFuncs_, invJ_(ii).arr, solidMesh4Sim.meshType);			
		Ke = ElementStiffMatrix(eleBmat, eleDmat, wgts, detJ_(:,ii));							
		eKs = Ke(eKk);
		sK(:,index) = eKs;						
	end						
	iK = eDofMat_(:,eKi)';
	jK = eDofMat_(:,eKj)';
	K_ = sparse(iK, jK, sK, numDOFs_, numDOFs_);
	K_ = K_ + K_' - diag(diag(K_));
	K_ = K_(freeDOFs_, freeDOFs_);
	
	%%Solving Equilibrium Equation
	climactericDOF = 5.0e5;
	if numDOFs_>climactericDOF
		deformation(freeDOFs_) = K_\Fvec;
	else
		[preconditionerL, preconditionerU] = ilu(K_);
		Preconditioning = @(x) preconditionerU\(preconditionerL\x);	
		ATX = @(x) K_*x;
		deformation(freeDOFs_) = myCG(ATX, Preconditioning, Fvec, 'printP_ON');	
	end
	
	%%Stress Analysis
	numNodsAroundEleVec_ = zeros(solidMesh4Sim.numNodes,1);
	for ii=1:solidMesh4Sim.numElements
		iEleNodes = solidMesh4Sim.eNodMat(ii,:);
		numNodsAroundEleVec_(iEleNodes,1) = numNodsAroundEleVec_(iEleNodes,1) + 1;
	end	
	cartesianStress = zeros(solidMesh4Sim.numNodes, 6);
	Ns = GetElementStressInterpolationMatrix(gaussIPs, solidMesh4Sim.meshType); OTP = inv(Ns);
	for ii=1:solidMesh4Sim.numElements
		iEleU = deformation(eDofMat_(ii,:),1);
		iMatrixB = ElementStrainMatrix(deShapeFuncs_, invJ_(ii).arr, solidMesh4Sim.meshType);
		stressGaussPoints = eleDmat * (iMatrixB*iEleU);
		stressNodes = OTP*stressGaussPoints;
		iNodes = solidMesh4Sim.eNodMat(ii,:);
		cartesianStress(iNodes,:) = reshape(stressNodes, 6, nEleNodes)' + cartesianStress(iNodes,:);
	end
	cartesianStress = cartesianStress./numNodsAroundEleVec_;
	
	compliance = deformation(freeDOFs_)' * (K_ * deformation(freeDOFs_));
end

function [eKi, eKj, eKk] = GetLowerEleStiffMatIndices(meshType)
	switch meshType
		case 'HEX', dimK = 24;
		case 'TET', dimK = 12;
	end
	rowMat = (1:dimK)'; rowMat = repmat(rowMat, 1, dimK);
	colMat = (1:dimK); colMat = repmat(colMat, dimK, 1);
	valMat = (1:dimK^2)'; valMat = reshape(valMat, dimK, dimK);
	eKi = tril(rowMat); [~, ~, eKi] = find(eKi);
	eKj = tril(colMat); [~, ~, eKj] = find(eKj);
	eKk = tril(valMat); [~, ~, eKk] = find(eKk);
end

function D = ElementElasticityMatrix(E, nu, meshType)	
	switch meshType
		case 'TET'
			HL = HookeLaw_SOLID(E, nu);
			D = zeros(24);
			for ii=1:4
				index = (ii-1)*6+1:ii*6;
				D(index,index) = HL;
			end			
		case 'HEX'
			HL = HookeLaw_SOLID(E, nu);
			D = zeros(48);
			for ii=1:8
				index = (ii-1)*6+1:ii*6;
				D(index,index) = HL;
			end						
	end
	D = sparse(D);
end

function HL = HookeLaw_SOLID(E, nu)
	cons1 = (1+nu)*(1-2*nu);
	cons2 = 2*(1+nu);
	HL = [
		E*(1-nu)/cons1	E*nu/cons1		E*nu/cons1		0			0		0
		E*nu/cons1		E*(1-nu)/cons1	E*nu/cons1		0			0		0
		E*nu/cons1		E*nu/cons1		E*(1-nu)/cons1	0			0		0
		0				0				0				E/cons2		0		0
		0				0				0				0			E/cons2	0
		0				0				0				0			0		E/cons2
	];
end

function val = GaussianIntegral(meshType)
	switch meshType
		case 'HEX'
			sqrt33 = sqrt(3)/3;
			val = [
				-sqrt33		sqrt33		sqrt33		-sqrt33		-sqrt33		sqrt33		sqrt33		-sqrt33
				-sqrt33		-sqrt33		sqrt33		sqrt33		-sqrt33		-sqrt33		sqrt33		sqrt33
				-sqrt33		-sqrt33		-sqrt33		-sqrt33		sqrt33		sqrt33		sqrt33		sqrt33
				1			1			1			1			1			1			1			1	
			];			
		case 'TET'
			alp = 0.58541020;
			bet = 0.13819660;	
			w = 1/6/4;		
			val = [
				bet		alp		bet		bet
				bet		bet		alp		bet
				bet		bet		bet		alp
				w		w		w		w
			];
	end
end

function dN = DeShapeFunction(paraCoords, meshType)	
	%% 	paraCoords = [
	%%		s1 s2 s3 ...
	%%		t1 t2 t3 ...
	%%		p1 p2 p3 ...
	%% ]
	switch meshType
		case 'TET'
			s = paraCoords(:,1);
			t = paraCoords(:,2);
			p = paraCoords(:,3);
			numVars = length(s); tmp = ones(numVars,1);
			
			dN1ds = -1*tmp; dN2ds = 1*tmp; dN3ds = 0*tmp; dN4ds = 0*tmp;
			dN1dt = -1*tmp; dN2dt = 0*tmp; dN3dt = 1*tmp; dN4dt = 0*tmp;	
			dN1dp = -1*tmp; dN2dp = 0*tmp; dN3dp = 0*tmp; dN4dp = 1*tmp;

			dN = zeros(3*length(s), 4);
			dN(1:3:end,:) = [dN1ds dN2ds dN3ds dN4ds];
			dN(2:3:end,:) = [dN1dt dN2dt dN3dt dN4dt];
			dN(3:3:end,:) = [dN1dp dN2dp dN3dp dN4dp];				
		case 'HEX'
			s = paraCoords(:,1);
			t = paraCoords(:,2);
			p = paraCoords(:,3);
			dN1ds = -0.125*(1-t).*(1-p); dN2ds = 0.125*(1-t).*(1-p); 
			dN3ds = 0.125*(1+t).*(1-p);  dN4ds = -0.125*(1+t).*(1-p);
			dN5ds = -0.125*(1-t).*(1+p); dN6ds = 0.125*(1-t).*(1+p); 
			dN7ds = 0.125*(1+t).*(1+p);  dN8ds = -0.125*(1+t).*(1+p);
		
			dN1dt = -0.125*(1-s).*(1-p); dN2dt = -0.125*(1+s).*(1-p); 
			dN3dt = 0.125*(1+s).*(1-p);  dN4dt = 0.125*(1-s).*(1-p);
			dN5dt = -0.125*(1-s).*(1+p); dN6dt = -0.125*(1+s).*(1+p); 
			dN7dt = 0.125*(1+s).*(1+p);  dN8dt = 0.125*(1-s).*(1+p);
		
			dN1dp = -0.125*(1-s).*(1-t); dN2dp = -0.125*(1+s).*(1-t); 
			dN3dp = -0.125*(1+s).*(1+t); dN4dp = -0.125*(1-s).*(1+t);
			dN5dp = 0.125*(1-s).*(1-t);  dN6dp = 0.125*(1+s).*(1-t); 
			dN7dp = 0.125*(1+s).*(1+t);  dN8dp = 0.125*(1-s).*(1+t);
			
			dN = zeros(3*length(s), 8);
			dN(1:3:end,:) = [dN1ds dN2ds dN3ds dN4ds dN5ds dN6ds dN7ds dN8ds];
			dN(2:3:end,:) = [dN1dt dN2dt dN3dt dN4dt dN5dt dN6dt dN7dt dN8dt];
			dN(3:3:end,:) = [dN1dp dN2dp dN3dp dN4dp dN5dp dN6dp dN7dp dN8dp];			
			
		
	end
end

function B = ElementStrainMatrix(dShape, invJ, meshType)
	derivatives = invJ * dShape;
	switch meshType
		case 'TET'
			dNds1 = derivatives(1,:);	dNdt1 = derivatives(2,:);	dNdp1 = derivatives(3,:);
			dNds2 = derivatives(4,:);	dNdt2 = derivatives(5,:);	dNdp2 = derivatives(6,:);
			dNds3 = derivatives(7,:);	dNdt3 = derivatives(8,:);	dNdp3 = derivatives(9,:);	
			dNds4 = derivatives(10,:);	dNdt4 = derivatives(11,:);	dNdp4 = derivatives(12,:);	
			
			B1_1 = [dNds1(1) 0 0; 0 dNdt1(1) 0; 0 0 dNdp1(1); 0 dNdp1(1) dNdt1(1); dNdp1(1) 0 dNds1(1); dNdt1(1) dNds1(1) 0];			
			B1_2 = [dNds1(2) 0 0; 0 dNdt1(2) 0; 0 0 dNdp1(2); 0 dNdp1(2) dNdt1(2); dNdp1(2) 0 dNds1(2); dNdt1(2) dNds1(2) 0];			
			B1_3 = [dNds1(3) 0 0; 0 dNdt1(3) 0; 0 0 dNdp1(3); 0 dNdp1(3) dNdt1(3); dNdp1(3) 0 dNds1(3); dNdt1(3) dNds1(3) 0];			
			B1_4 = [dNds1(4) 0 0; 0 dNdt1(4) 0; 0 0 dNdp1(4); 0 dNdp1(4) dNdt1(4); dNdp1(4) 0 dNds1(4); dNdt1(4) dNds1(4) 0];		
		
			B2_1 = [dNds2(1) 0 0; 0 dNdt2(1) 0; 0 0 dNdp2(1); 0 dNdp2(1) dNdt2(1); dNdp2(1) 0 dNds2(1); dNdt2(1) dNds2(1) 0];			
			B2_2 = [dNds2(2) 0 0; 0 dNdt2(2) 0; 0 0 dNdp2(2); 0 dNdp2(2) dNdt2(2); dNdp2(2) 0 dNds2(2); dNdt2(2) dNds2(2) 0];			
			B2_3 = [dNds2(3) 0 0; 0 dNdt2(3) 0; 0 0 dNdp2(3); 0 dNdp2(3) dNdt2(3); dNdp2(3) 0 dNds2(3); dNdt2(3) dNds2(3) 0];			
			B2_4 = [dNds2(4) 0 0; 0 dNdt2(4) 0; 0 0 dNdp2(4); 0 dNdp2(4) dNdt2(4); dNdp2(4) 0 dNds2(4); dNdt2(4) dNds2(4) 0];			
		
			B3_1 = [dNds3(1) 0 0; 0 dNdt3(1) 0; 0 0 dNdp3(1); 0 dNdp3(1) dNdt3(1); dNdp3(1) 0 dNds3(1); dNdt3(1) dNds3(1) 0];			
			B3_2 = [dNds3(2) 0 0; 0 dNdt3(2) 0; 0 0 dNdp3(2); 0 dNdp3(2) dNdt3(2); dNdp3(2) 0 dNds3(2); dNdt3(2) dNds3(2) 0];			
			B3_3 = [dNds3(3) 0 0; 0 dNdt3(3) 0; 0 0 dNdp3(3); 0 dNdp3(3) dNdt3(3); dNdp3(3) 0 dNds3(3); dNdt3(3) dNds3(3) 0];			
			B3_4 = [dNds3(4) 0 0; 0 dNdt3(4) 0; 0 0 dNdp3(4); 0 dNdp3(4) dNdt3(4); dNdp3(4) 0 dNds3(4); dNdt3(4) dNds3(4) 0];
		
			B4_1 = [dNds4(1) 0 0; 0 dNdt4(1) 0; 0 0 dNdp4(1); 0 dNdp4(1) dNdt4(1); dNdp4(1) 0 dNds4(1); dNdt4(1) dNds4(1) 0];			
			B4_2 = [dNds4(2) 0 0; 0 dNdt4(2) 0; 0 0 dNdp4(2); 0 dNdp4(2) dNdt4(2); dNdp4(2) 0 dNds4(2); dNdt4(2) dNds4(2) 0];			
			B4_3 = [dNds4(3) 0 0; 0 dNdt4(3) 0; 0 0 dNdp4(3); 0 dNdp4(3) dNdt4(3); dNdp4(3) 0 dNds4(3); dNdt4(3) dNds4(3) 0];			                                                                                                      
			B4_4 = [dNds4(4) 0 0; 0 dNdt4(4) 0; 0 0 dNdp4(4); 0 dNdp4(4) dNdt4(4); dNdp4(4) 0 dNds4(4); dNdt4(4) dNds4(4) 0];

			B = [
				B1_1 B1_2 B1_3 B1_4
				B2_1 B2_2 B2_3 B2_4
				B3_1 B3_2 B3_3 B3_4
				B4_1 B4_2 B4_3 B4_4	
			];
		case 'HEX'
			dNds1 = derivatives(1,:);	dNdt1 = derivatives(2,:);	dNdp1 = derivatives(3,:);
			dNds2 = derivatives(4,:);	dNdt2 = derivatives(5,:);	dNdp2 = derivatives(6,:);
			dNds3 = derivatives(7,:);	dNdt3 = derivatives(8,:);	dNdp3 = derivatives(9,:);	
			dNds4 = derivatives(10,:);	dNdt4 = derivatives(11,:);	dNdp4 = derivatives(12,:);	
			dNds5 = derivatives(13,:);	dNdt5 = derivatives(14,:);	dNdp5 = derivatives(15,:);	
			dNds6 = derivatives(16,:);	dNdt6 = derivatives(17,:);	dNdp6 = derivatives(18,:);	
			dNds7 = derivatives(19,:);	dNdt7 = derivatives(20,:);	dNdp7 = derivatives(21,:);
			dNds8 = derivatives(22,:);	dNdt8 = derivatives(23,:);	dNdp8 = derivatives(24,:);
			
			B1_1 = [dNds1(1) 0 0; 0 dNdt1(1) 0; 0 0 dNdp1(1); 0 dNdp1(1) dNdt1(1); dNdp1(1) 0 dNds1(1); dNdt1(1) dNds1(1) 0];			
			B1_2 = [dNds1(2) 0 0; 0 dNdt1(2) 0; 0 0 dNdp1(2); 0 dNdp1(2) dNdt1(2); dNdp1(2) 0 dNds1(2); dNdt1(2) dNds1(2) 0];			
			B1_3 = [dNds1(3) 0 0; 0 dNdt1(3) 0; 0 0 dNdp1(3); 0 dNdp1(3) dNdt1(3); dNdp1(3) 0 dNds1(3); dNdt1(3) dNds1(3) 0];			
			B1_4 = [dNds1(4) 0 0; 0 dNdt1(4) 0; 0 0 dNdp1(4); 0 dNdp1(4) dNdt1(4); dNdp1(4) 0 dNds1(4); dNdt1(4) dNds1(4) 0];			
			B1_5 = [dNds1(5) 0 0; 0 dNdt1(5) 0; 0 0 dNdp1(5); 0 dNdp1(5) dNdt1(5); dNdp1(5) 0 dNds1(5); dNdt1(5) dNds1(5) 0];			
			B1_6 = [dNds1(6) 0 0; 0 dNdt1(6) 0; 0 0 dNdp1(6); 0 dNdp1(6) dNdt1(6); dNdp1(6) 0 dNds1(6); dNdt1(6) dNds1(6) 0];			
			B1_7 = [dNds1(7) 0 0; 0 dNdt1(7) 0; 0 0 dNdp1(7); 0 dNdp1(7) dNdt1(7); dNdp1(7) 0 dNds1(7); dNdt1(7) dNds1(7) 0];			
			B1_8 = [dNds1(8) 0 0; 0 dNdt1(8) 0; 0 0 dNdp1(8); 0 dNdp1(8) dNdt1(8); dNdp1(8) 0 dNds1(8); dNdt1(8) dNds1(8) 0];
		
			B2_1 = [dNds2(1) 0 0; 0 dNdt2(1) 0; 0 0 dNdp2(1); 0 dNdp2(1) dNdt2(1); dNdp2(1) 0 dNds2(1); dNdt2(1) dNds2(1) 0];			
			B2_2 = [dNds2(2) 0 0; 0 dNdt2(2) 0; 0 0 dNdp2(2); 0 dNdp2(2) dNdt2(2); dNdp2(2) 0 dNds2(2); dNdt2(2) dNds2(2) 0];			
			B2_3 = [dNds2(3) 0 0; 0 dNdt2(3) 0; 0 0 dNdp2(3); 0 dNdp2(3) dNdt2(3); dNdp2(3) 0 dNds2(3); dNdt2(3) dNds2(3) 0];			
			B2_4 = [dNds2(4) 0 0; 0 dNdt2(4) 0; 0 0 dNdp2(4); 0 dNdp2(4) dNdt2(4); dNdp2(4) 0 dNds2(4); dNdt2(4) dNds2(4) 0];			
			B2_5 = [dNds2(5) 0 0; 0 dNdt2(5) 0; 0 0 dNdp2(5); 0 dNdp2(5) dNdt2(5); dNdp2(5) 0 dNds2(5); dNdt2(5) dNds2(5) 0];			
			B2_6 = [dNds2(6) 0 0; 0 dNdt2(6) 0; 0 0 dNdp2(6); 0 dNdp2(6) dNdt2(6); dNdp2(6) 0 dNds2(6); dNdt2(6) dNds2(6) 0];			
			B2_7 = [dNds2(7) 0 0; 0 dNdt2(7) 0; 0 0 dNdp2(7); 0 dNdp2(7) dNdt2(7); dNdp2(7) 0 dNds2(7); dNdt2(7) dNds2(7) 0];			
			B2_8 = [dNds2(8) 0 0; 0 dNdt2(8) 0; 0 0 dNdp2(8); 0 dNdp2(8) dNdt2(8); dNdp2(8) 0 dNds2(8); dNdt2(8) dNds2(8) 0];			
		
			B3_1 = [dNds3(1) 0 0; 0 dNdt3(1) 0; 0 0 dNdp3(1); 0 dNdp3(1) dNdt3(1); dNdp3(1) 0 dNds3(1); dNdt3(1) dNds3(1) 0];			
			B3_2 = [dNds3(2) 0 0; 0 dNdt3(2) 0; 0 0 dNdp3(2); 0 dNdp3(2) dNdt3(2); dNdp3(2) 0 dNds3(2); dNdt3(2) dNds3(2) 0];			
			B3_3 = [dNds3(3) 0 0; 0 dNdt3(3) 0; 0 0 dNdp3(3); 0 dNdp3(3) dNdt3(3); dNdp3(3) 0 dNds3(3); dNdt3(3) dNds3(3) 0];			
			B3_4 = [dNds3(4) 0 0; 0 dNdt3(4) 0; 0 0 dNdp3(4); 0 dNdp3(4) dNdt3(4); dNdp3(4) 0 dNds3(4); dNdt3(4) dNds3(4) 0];			
			B3_5 = [dNds3(5) 0 0; 0 dNdt3(5) 0; 0 0 dNdp3(5); 0 dNdp3(5) dNdt3(5); dNdp3(5) 0 dNds3(5); dNdt3(5) dNds3(5) 0];			
			B3_6 = [dNds3(6) 0 0; 0 dNdt3(6) 0; 0 0 dNdp3(6); 0 dNdp3(6) dNdt3(6); dNdp3(6) 0 dNds3(6); dNdt3(6) dNds3(6) 0];			
			B3_7 = [dNds3(7) 0 0; 0 dNdt3(7) 0; 0 0 dNdp3(7); 0 dNdp3(7) dNdt3(7); dNdp3(7) 0 dNds3(7); dNdt3(7) dNds3(7) 0];			
			B3_8 = [dNds3(8) 0 0; 0 dNdt3(8) 0; 0 0 dNdp3(8); 0 dNdp3(8) dNdt3(8); dNdp3(8) 0 dNds3(8); dNdt3(8) dNds3(8) 0];
		
			B4_1 = [dNds4(1) 0 0; 0 dNdt4(1) 0; 0 0 dNdp4(1); 0 dNdp4(1) dNdt4(1); dNdp4(1) 0 dNds4(1); dNdt4(1) dNds4(1) 0];			
			B4_2 = [dNds4(2) 0 0; 0 dNdt4(2) 0; 0 0 dNdp4(2); 0 dNdp4(2) dNdt4(2); dNdp4(2) 0 dNds4(2); dNdt4(2) dNds4(2) 0];			
			B4_3 = [dNds4(3) 0 0; 0 dNdt4(3) 0; 0 0 dNdp4(3); 0 dNdp4(3) dNdt4(3); dNdp4(3) 0 dNds4(3); dNdt4(3) dNds4(3) 0];			                                                                                                      
			B4_4 = [dNds4(4) 0 0; 0 dNdt4(4) 0; 0 0 dNdp4(4); 0 dNdp4(4) dNdt4(4); dNdp4(4) 0 dNds4(4); dNdt4(4) dNds4(4) 0];			                                                                                                      
			B4_5 = [dNds4(5) 0 0; 0 dNdt4(5) 0; 0 0 dNdp4(5); 0 dNdp4(5) dNdt4(5); dNdp4(5) 0 dNds4(5); dNdt4(5) dNds4(5) 0];			                                                                                                      
			B4_6 = [dNds4(6) 0 0; 0 dNdt4(6) 0; 0 0 dNdp4(6); 0 dNdp4(6) dNdt4(6); dNdp4(6) 0 dNds4(6); dNdt4(6) dNds4(6) 0];					  		                                                                                      
			B4_7 = [dNds4(7) 0 0; 0 dNdt4(7) 0; 0 0 dNdp4(7); 0 dNdp4(7) dNdt4(7); dNdp4(7) 0 dNds4(7); dNdt4(7) dNds4(7) 0];					  		                                                                                      
			B4_8 = [dNds4(8) 0 0; 0 dNdt4(8) 0; 0 0 dNdp4(8); 0 dNdp4(8) dNdt4(8); dNdp4(8) 0 dNds4(8); dNdt4(8) dNds4(8) 0];
		
			B5_1 = [dNds5(1) 0 0; 0 dNdt5(1) 0; 0 0 dNdp5(1); 0 dNdp5(1) dNdt5(1); dNdp5(1) 0 dNds5(1); dNdt5(1) dNds5(1) 0];					  		 		 	                                                                          
			B5_2 = [dNds5(2) 0 0; 0 dNdt5(2) 0; 0 0 dNdp5(2); 0 dNdp5(2) dNdt5(2); dNdp5(2) 0 dNds5(2); dNdt5(2) dNds5(2) 0];					  		 		 	                                                                          
			B5_3 = [dNds5(3) 0 0; 0 dNdt5(3) 0; 0 0 dNdp5(3); 0 dNdp5(3) dNdt5(3); dNdp5(3) 0 dNds5(3); dNdt5(3) dNds5(3) 0];					  		 		 	                                                                          
			B5_4 = [dNds5(4) 0 0; 0 dNdt5(4) 0; 0 0 dNdp5(4); 0 dNdp5(4) dNdt5(4); dNdp5(4) 0 dNds5(4); dNdt5(4) dNds5(4) 0];					  		 		 	                                                                          
			B5_5 = [dNds5(5) 0 0; 0 dNdt5(5) 0; 0 0 dNdp5(5); 0 dNdp5(5) dNdt5(5); dNdp5(5) 0 dNds5(5); dNdt5(5) dNds5(5) 0];					  		 		 	                                                                          
			B5_6 = [dNds5(6) 0 0; 0 dNdt5(6) 0; 0 0 dNdp5(6); 0 dNdp5(6) dNdt5(6); dNdp5(6) 0 dNds5(6); dNdt5(6) dNds5(6) 0];					  		 		 	                                                                          
			B5_7 = [dNds5(7) 0 0; 0 dNdt5(7) 0; 0 0 dNdp5(7); 0 dNdp5(7) dNdt5(7); dNdp5(7) 0 dNds5(7); dNdt5(7) dNds5(7) 0];					  		 		 	                                                                          
			B5_8 = [dNds5(8) 0 0; 0 dNdt5(8) 0; 0 0 dNdp5(8); 0 dNdp5(8) dNdt5(8); dNdp5(8) 0 dNds5(8); dNdt5(8) dNds5(8) 0];
		
			B6_1 = [dNds6(1) 0 0; 0 dNdt6(1) 0; 0 0 dNdp6(1); 0 dNdp6(1) dNdt6(1); dNdp6(1) 0 dNds6(1); dNdt6(1) dNds6(1) 0];					  		 		 	 	                                                                      
			B6_2 = [dNds6(2) 0 0; 0 dNdt6(2) 0; 0 0 dNdp6(2); 0 dNdp6(2) dNdt6(2); dNdp6(2) 0 dNds6(2); dNdt6(2) dNds6(2) 0];					  		 		 	 	 	                                                                  
			B6_3 = [dNds6(3) 0 0; 0 dNdt6(3) 0; 0 0 dNdp6(3); 0 dNdp6(3) dNdt6(3); dNdp6(3) 0 dNds6(3); dNdt6(3) dNds6(3) 0];					  		 		 	 	 			                                                          
			B6_4 = [dNds6(4) 0 0; 0 dNdt6(4) 0; 0 0 dNdp6(4); 0 dNdp6(4) dNdt6(4); dNdp6(4) 0 dNds6(4); dNdt6(4) dNds6(4) 0];					  		 		 	 	 		                                                              
			B6_5 = [dNds6(5) 0 0; 0 dNdt6(5) 0; 0 0 dNdp6(5); 0 dNdp6(5) dNdt6(5); dNdp6(5) 0 dNds6(5); dNdt6(5) dNds6(5) 0];					  		 		 	 	 			                                                          
			B6_6 = [dNds6(6) 0 0; 0 dNdt6(6) 0; 0 0 dNdp6(6); 0 dNdp6(6) dNdt6(6); dNdp6(6) 0 dNds6(6); dNdt6(6) dNds6(6) 0];					  		 		 	 	 			                                                          
			B6_7 = [dNds6(7) 0 0; 0 dNdt6(7) 0; 0 0 dNdp6(7); 0 dNdp6(7) dNdt6(7); dNdp6(7) 0 dNds6(7); dNdt6(7) dNds6(7) 0];					  		 		 	 	 			                                                          
			B6_8 = [dNds6(8) 0 0; 0 dNdt6(8) 0; 0 0 dNdp6(8); 0 dNdp6(8) dNdt6(8); dNdp6(8) 0 dNds6(8); dNdt6(8) dNds6(8) 0];					  		 		 	 	 			                                                          
			
			B7_1 = [dNds7(1) 0 0; 0 dNdt7(1) 0; 0 0 dNdp7(1); 0 dNdp7(1) dNdt7(1); dNdp7(1) 0 dNds7(1); dNdt7(1) dNds7(1) 0];					  		 		 	 	 			                                                          
			B7_2 = [dNds7(2) 0 0; 0 dNdt7(2) 0; 0 0 dNdp7(2); 0 dNdp7(2) dNdt7(2); dNdp7(2) 0 dNds7(2); dNdt7(2) dNds7(2) 0];					  		 		 	 	 			                                                          
			B7_3 = [dNds7(3) 0 0; 0 dNdt7(3) 0; 0 0 dNdp7(3); 0 dNdp7(3) dNdt7(3); dNdp7(3) 0 dNds7(3); dNdt7(3) dNds7(3) 0];					  		 		 	 	 			                                                          
			B7_4 = [dNds7(4) 0 0; 0 dNdt7(4) 0; 0 0 dNdp7(4); 0 dNdp7(4) dNdt7(4); dNdp7(4) 0 dNds7(4); dNdt7(4) dNds7(4) 0];					  		 		 	 	 			                                                          
			B7_5 = [dNds7(5) 0 0; 0 dNdt7(5) 0; 0 0 dNdp7(5); 0 dNdp7(5) dNdt7(5); dNdp7(5) 0 dNds7(5); dNdt7(5) dNds7(5) 0];					  		 		 	 	 			                                                            
			B7_6 = [dNds7(6) 0 0; 0 dNdt7(6) 0; 0 0 dNdp7(6); 0 dNdp7(6) dNdt7(6); dNdp7(6) 0 dNds7(6); dNdt7(6) dNds7(6) 0];					  		 		 	 	 			                                                            
			B7_7 = [dNds7(7) 0 0; 0 dNdt7(7) 0; 0 0 dNdp7(7); 0 dNdp7(7) dNdt7(7); dNdp7(7) 0 dNds7(7); dNdt7(7) dNds7(7) 0];					  		 		 	 	 			                                                            
			B7_8 = [dNds7(8) 0 0; 0 dNdt7(8) 0; 0 0 dNdp7(8); 0 dNdp7(8) dNdt7(8); dNdp7(8) 0 dNds7(8); dNdt7(8) dNds7(8) 0];
			
			B8_1 = [dNds8(1) 0 0; 0 dNdt8(1) 0; 0 0 dNdp8(1); 0 dNdp8(1) dNdt8(1); dNdp8(1) 0 dNds8(1); dNdt8(1) dNds8(1) 0];					  		 		 	 	 			                                                            
			B8_2 = [dNds8(2) 0 0; 0 dNdt8(2) 0; 0 0 dNdp8(2); 0 dNdp8(2) dNdt8(2); dNdp8(2) 0 dNds8(2); dNdt8(2) dNds8(2) 0];					  		 		 	 	 			                                                            
			B8_3 = [dNds8(3) 0 0; 0 dNdt8(3) 0; 0 0 dNdp8(3); 0 dNdp8(3) dNdt8(3); dNdp8(3) 0 dNds8(3); dNdt8(3) dNds8(3) 0];					  		 		 	 	 			                                                            
			B8_4 = [dNds8(4) 0 0; 0 dNdt8(4) 0; 0 0 dNdp8(4); 0 dNdp8(4) dNdt8(4); dNdp8(4) 0 dNds8(4); dNdt8(4) dNds8(4) 0];					  		 		 	 	 			                                                            
			B8_5 = [dNds8(5) 0 0; 0 dNdt8(5) 0; 0 0 dNdp8(5); 0 dNdp8(5) dNdt8(5); dNdp8(5) 0 dNds8(5); dNdt8(5) dNds8(5) 0];					  		 		 	 	 			                                                            
			B8_6 = [dNds8(6) 0 0; 0 dNdt8(6) 0; 0 0 dNdp8(6); 0 dNdp8(6) dNdt8(6); dNdp8(6) 0 dNds8(6); dNdt8(6) dNds8(6) 0];					  		 		 	 	 			                                                            
			B8_7 = [dNds8(7) 0 0; 0 dNdt8(7) 0; 0 0 dNdp8(7); 0 dNdp8(7) dNdt8(7); dNdp8(7) 0 dNds8(7); dNdt8(7) dNds8(7) 0];					  		 		 	 	 			                                                            
			B8_8 = [dNds8(8) 0 0; 0 dNdt8(8) 0; 0 0 dNdp8(8); 0 dNdp8(8) dNdt8(8); dNdp8(8) 0 dNds8(8); dNdt8(8) dNds8(8) 0];
		
			B = [
				B1_1 B1_2 B1_3 B1_4 B1_5 B1_6 B1_7 B1_8
				B2_1 B2_2 B2_3 B2_4 B2_5 B2_6 B2_7 B2_8
				B3_1 B3_2 B3_3 B3_4 B3_5 B3_6 B3_7 B3_8
				B4_1 B4_2 B4_3 B4_4 B4_5 B4_6 B4_7 B4_8
				B5_1 B5_2 B5_3 B5_4 B5_5 B5_6 B5_7 B5_8
				B6_1 B6_2 B6_3 B6_4 B6_5 B6_6 B6_7 B6_8
				B7_1 B7_2 B7_3 B7_4 B7_5 B7_6 B7_7 B7_8
				B8_1 B8_2 B8_3 B8_4 B8_5 B8_6 B8_7 B8_8		
			];		
	end
end

function Ke = ElementStiffMatrix(B, D, w, detJ)
	wgt = w.*detJ;	wgt = repmat(wgt, 1, 6);
	wgt = reshape(wgt', 1, numel(wgt));
	Ke = B'*(D.*wgt)*B;
	%%ref
	%Ke = B'*D*sparse(diag(wgt))*B;	
end

function x = myCG(ATX, Preconditioning, b, printP, varargin)
	tol = 1.0e-3; maxIT = 20000;
	n = length(b);
	normB = norm(b);
	its = 0;
	if 5==nargin
		x = varargin{1};
	else
		x = zeros(n,1);
	end
	r1 = b - ATX(x);
	z1 = zeros(n,1);

	while its <= maxIT
		its = its + 1;
		z2 = Preconditioning(r1);
		if 1==its
			p2 = z2;
		else
			beta = r1'*z2/(r0'*z1);
			p2 = z2 + beta*p1;			
		end
		valMTV = ATX(p2);	
		alpha = r1'*z2/(p2'*valMTV);
		x = x + alpha*p2;		
		r2 = r1 - alpha*valMTV;
		
		resnorm = norm(r2)/normB;
		if strcmp(printP, 'printP_ON')
			disp([' It.: ' sprintf('%4i',its) ' Res.: ' sprintf('%16.6e',resnorm)]);
		end
		if resnorm<tol
			disp(['Conjugate Gradient Solver Converged at Iteration' sprintf('%5i', its) ' to a Solution with Relative Residual' ...
					sprintf('%16.6e',resnorm)]);
			break;
		end		
		%%update
		z1 = z2;
		p1 = p2;
		r0 = r1;
		r1 = r2;
	end	

	if its > maxIT
		warning('Exceed the maximum iterate numbers');
		disp(['The iterative process stops at residual = ' sprintf('%10.4f',norm(r2))]);		
	end
end

function Ns = GetElementStressInterpolationMatrix(gaussIPs, meshType)
	switch meshType 			
		case 'TET'
			N = ShapeFunction(gaussIPs, meshType);
			Ns = sparse(24,24);
			ii = 6*(1:4);
			Ns(1,ii-5) = N(1,:); Ns(2,ii-4) = N(1,:); Ns(3,ii-3) = N(1,:);
			Ns(4,ii-2) = N(1,:); Ns(5,ii-1) = N(1,:); Ns(6,ii) = N(1,:);
			
			Ns(7,ii-5) = N(2,:); Ns(8,ii-4) = N(2,:); Ns(9,ii-3) = N(2,:);
			Ns(10,ii-2) = N(2,:); Ns(11,ii-1) = N(2,:); Ns(12,ii) = N(2,:);

			Ns(13,ii-5) = N(3,:); Ns(14,ii-4) = N(3,:); Ns(15,ii-3) = N(3,:);
			Ns(16,ii-2) = N(3,:); Ns(17,ii-1) = N(3,:); Ns(18,ii) = N(3,:);	

			Ns(19,ii-5) = N(4,:); Ns(20,ii-4) = N(4,:); Ns(21,ii-3) = N(4,:);
			Ns(22,ii-2) = N(4,:); Ns(23,ii-1) = N(4,:); Ns(24,ii) = N(4,:);					
		case 'HEX'
			N = ShapeFunction(gaussIPs, meshType);
			Ns = sparse(48,48);
			ii = 6*(1:8);
			Ns(1,ii-5) = N(1,:); Ns(2,ii-4) = N(1,:); Ns(3,ii-3) = N(1,:);
			Ns(4,ii-2) = N(1,:); Ns(5,ii-1) = N(1,:); Ns(6,ii) = N(1,:);
			
			Ns(7,ii-5) = N(2,:); Ns(8,ii-4) = N(2,:); Ns(9,ii-3) = N(2,:);
			Ns(10,ii-2) = N(2,:); Ns(11,ii-1) = N(2,:); Ns(12,ii) = N(2,:);

			Ns(13,ii-5) = N(3,:); Ns(14,ii-4) = N(3,:); Ns(15,ii-3) = N(3,:);
			Ns(16,ii-2) = N(3,:); Ns(17,ii-1) = N(3,:); Ns(18,ii) = N(3,:);	

			Ns(19,ii-5) = N(4,:); Ns(20,ii-4) = N(4,:); Ns(21,ii-3) = N(4,:);
			Ns(22,ii-2) = N(4,:); Ns(23,ii-1) = N(4,:); Ns(24,ii) = N(4,:);

			Ns(25,ii-5) = N(5,:); Ns(26,ii-4) = N(5,:); Ns(27,ii-3) = N(5,:);
			Ns(28,ii-2) = N(5,:); Ns(29,ii-1) = N(5,:); Ns(30,ii) = N(5,:);	

			Ns(31,ii-5) = N(6,:); Ns(32,ii-4) = N(6,:); Ns(33,ii-3) = N(6,:);
			Ns(34,ii-2) = N(6,:); Ns(35,ii-1) = N(6,:); Ns(36,ii) = N(6,:);	

			Ns(37,ii-5) = N(7,:); Ns(38,ii-4) = N(7,:); Ns(39,ii-3) = N(7,:);
			Ns(40,ii-2) = N(7,:); Ns(41,ii-1) = N(7,:); Ns(42,ii) = N(7,:);	

			Ns(43,ii-5) = N(8,:); Ns(44,ii-4) = N(8,:); Ns(45,ii-3) = N(8,:);
			Ns(46,ii-2) = N(8,:); Ns(47,ii-1) = N(8,:); Ns(48,ii) = N(8,:);				
	end
end

function N = ShapeFunction(paraCoords, meshType)
	%% 	paraCoords = [
	%%		s1 s2 s3 ...
	%%		t1 t2 t3 ...
	%%		p1 p2 p3 ...
	%% ]
	switch meshType	
		case 'TET'
			s = paraCoords(:,1);
			t = paraCoords(:,2);
			p = paraCoords(:,3);
			N = zeros(length(s), 4);
			N(:,1) = 1-s-t-p;
			N(:,2) = s;
			N(:,3) = t;
			N(:,4) = p;
		case 'HEX'
			%				*8			*7
			%			*5			*6
			%					p
			%				   |__s (parametric coordinate system)
			%				  /-t
			%				*4			*3
			%			*1			*2
			%
			%				nodes
			s = paraCoords(:,1);
			t = paraCoords(:,2);
			p = paraCoords(:,3);
			N = zeros(length(s), 8);
			N(:,1) = 0.125*(1-s).*(1-t).*(1-p);
			N(:,2) = 0.125*(1+s).*(1-t).*(1-p);
			N(:,3) = 0.125*(1+s).*(1+t).*(1-p);
			N(:,4) = 0.125*(1-s).*(1+t).*(1-p);
			N(:,5) = 0.125*(1-s).*(1-t).*(1+p);
			N(:,6) = 0.125*(1+s).*(1-t).*(1+p);
			N(:,7) = 0.125*(1+s).*(1+t).*(1+p);
			N(:,8) = 0.125*(1-s).*(1+t).*(1+p);						
	end	
end
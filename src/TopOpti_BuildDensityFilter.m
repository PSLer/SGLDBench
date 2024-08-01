function TopOpti_BuildDensityFilter()
	global outPath_;
	global meshHierarchy_;
	global rMin_;
	global H_;
	global Hs_;
	
	%%1. setup
	rMin = rMin_;
	eleSize = meshHierarchy_(1).eleSize(1);
	numElements = meshHierarchy_(1).numElements;
	eleIndexMapForward = meshHierarchy_(1).eleMapForward;	
	p = 3;
	eleCentroidList = niftiread(strcat(outPath_, 'cache_eleCentroidList.nii')); eleCentroidList = double(eleCentroidList);
	
	%%2. Relate Element Adjacency
	iH = zeros(numElements*(2*(ceil(rMin)-1)+1)^p, 1, 'int32');
	jH = iH;
	iIndex = 0;

	%%	1	4	7		10	 13	  16		19	 22	  25
	%%	2	5	8		11	 14*  17		20	 23	  26
	%%	3	6	9		12	 15   18		21	 24	  27
	%%	 bottom				middle				top
	resX = meshHierarchy_(1).resX;
	resY = meshHierarchy_(1).resY;
	resZ = meshHierarchy_(1).resZ;
	
	for kk = 1:resZ
		for ii = 1:resX
			for jj = 1:resY			
				e1 = eleIndexMapForward((kk-1)*resX*resY+(ii-1)*resY+jj);
				if e1
					for kk2 = max(kk-(ceil(rMin)-1),1):min(kk+(ceil(rMin)-1),resZ)
						for ii2 = max(ii-(ceil(rMin)-1),1):min(ii+(ceil(rMin)-1),resX)
							for jj2 = max(jj-(ceil(rMin)-1),1):min(jj+(ceil(rMin)-1),resY)
								e2 = eleIndexMapForward((kk2-1)*resX*resY+(ii2-1)*resY+jj2);
								if e2
									iIndex = iIndex + 1;
									iH(iIndex) = e1;
									jH(iIndex) = e2;									
								end
							end
						end
					end
				end
			end
		end
	end	

	iH = iH(1:iIndex,:);
	jH = jH(1:iIndex,:);
	
	%%3. Construct Density Filter
	if 0 %%minimum thickness control
		rr = 2*vecnorm(eleCentroidList(iH,:)-eleCentroidList(jH,:),2,2)/(rMin*eleSize);
		rr2 = rr.^2;
		sH = 1-6*rr2+8*rr.*rr2-3*rr2.^2;			
	else %%remove checkerboard pattern
		sH = rMin*eleSize - vecnorm(eleCentroidList(iH,:)-eleCentroidList(jH,:),2,2);
		sH(sH<0) = 0;		
	end

	blockIndex = Solving_MissionPartition(size(iH,1), 1.0e8);
	H_ = sparse(numElements, numElements);
	for ii=1:size(blockIndex,1)
		rangeIndex = (blockIndex(ii,1):blockIndex(ii,2))';
		tmpH = sparse(iH(rangeIndex,:), jH(rangeIndex,:), sH(rangeIndex,:), numElements, numElements);					
		H_ = H_ + tmpH;
	end	
	Hs_ = sum(H_,2);	
end
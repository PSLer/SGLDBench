function Hs = TopOpti_SetupDensityFiltering(rMin, resolutions, numElements, eleMapForward)
	resX = resolutions(1); resY = resolutions(2); resZ = resolutions(3);
	Hs = zeros(numElements,1);
	for kk = 1:resZ
		for ii = 1:resX
			for jj = 1:resY			
				e1MapBack = (kk-1)*resX*resY+(ii-1)*resY+jj;
				e1 = eleMapForward(e1MapBack);
				if e1
					weightsSum = 0;
					for kk2 = max(kk-(ceil(rMin)-1),1):min(kk+(ceil(rMin)-1),resZ)
						for ii2 = max(ii-(ceil(rMin)-1),1):min(ii+(ceil(rMin)-1),resX)
							for jj2 = max(jj-(ceil(rMin)-1),1):min(jj+(ceil(rMin)-1),resY)
								e2MapBack = (kk2-1)*resX*resY+(ii2-1)*resY+jj2;
								e2 = eleMapForward(e2MapBack);
								if e2
									e2e1Weight = max(0.0, rMin - sqrt((ii-ii2)^2 + (jj-jj2)^2 + (kk-kk2)^2));
									% e1Weights(iIndex,1) = max(0.0, sqrt((ii-ii2)^2 + (jj-jj2)^2 + (kk-kk2)^2));
									% e1Val = e1Val + xList(e2) * e2e1Weight;
									weightsSum = weightsSum + e2e1Weight;
								end
							end
						end
					end
					Hs(e1) = weightsSum;				
				end
			end
		end
	end	
end
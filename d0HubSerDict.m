function dHub = d0HubSerDict(allP)
%DHUBCLASSREPDICT form Hub d map - priors over HI, SI so that the 
%  alpha and beta are each at least 1 and pH0, pS0 are the modes.
%  allP must include, for use here: pH0, pS0, dInitEv, initEvRat, resNHub, noiseFloor

pH0 = allP.pH0; pS0 = allP.pS0;   
try   % early versions of the code did not have initEvRat , the dInit H/S ratio.
   r = allP.initEvRat; 
   dInitH = allP.dInitEv * r / (1+r);
   dInitS = allP.dInitEv / (1+r); 
catch
   warning('about to use allP.initEv for both dInitH and S'); 
   dInitH = allP.initEv;  
   dInitS = dInitH;
end
resNHub=allP.resNHub;   
noiseFloor = allP.noiseFloor;

% dHub{1} for HI: 
alphaHI = pH0 * dInitH; 
betaHI = dInitH - alphaHI;  
cbd = [0 betacdf( (1:(resNHub-1))/resNHub, alphaHI, betaHI ) 1];
pmf = cbd(2:end) - cbd(1:resNHub) + noiseFloor; 
dHub{1} = (pmf' / sum(pmf)) * dInitH;
  
% dHub{2} for SI: 
alphaSI = pS0 * dInitS; 
betaSI = dInitS - alphaSI;  
cbd = [0 betacdf( (1:(resNHub-1))/resNHub, alphaSI, betaSI ) 1];
pmf = cbd(2:end) - cbd(1:resNHub) + noiseFloor; 
dHub{2} = (pmf' / sum(pmf)) * dInitS;

return;


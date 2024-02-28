function dHub = dHubClassSerDict(allP)
%DHUBCLASSREPDICT form Hub d map - priors over HI, SI so that the 
%  alpha and beta are each at least 1 and pH0, pS0 are the modes.
%  allP must include, for use here: pH0, pS0, dInitEv, initEvRat, resNHub, noiseFloor

if length(allP.HIv0) ~= 2
    error('allP.HIv0) ~= 2 not catered for yet');
end

pH0 = allP.pH0; pS0 = allP.pS0;   
try Upop = allP.Upop;  catch, Upop = 1;  end
try   % early versions of the code did not have initEvRat , the dInit S/H ratio.
   r = allP.initEvRat; 
   dInitH = allP.dInitEv * r / (1+r);
   dInitS = allP.dInitEv / (1+r); 
catch
   warning('about to use allP.initEv'); 
   dInitH = allP.initEv;  
   dInitS = dInitH;
end
% resNHub=allP.resNHub;   
pLeast = allP.noiseFloor;

frH = (pH0 - allP.HIv0(1)) / (allP.HIv0(2) - allP.HIv0(1));
if frH <= pLeast 
    frH = pLeast; 
elseif frH >= (1 - pLeast)
    frH = 1 - pLeast;
end
frS = (pS0 - allP.SIv0(1)) / (allP.SIv0(2) - allP.SIv0(1));
if frS <= pLeast 
    frS = pLeast; 
elseif frS >= (1-pLeast)
    frS = 1-pLeast;
end

dHub{1} = [ (1-frH)*(1+dInitH); frH * (1+dInitH)];
dHub{2} = [ (1-frS)*(1+dInitS); frS * (1+dInitS)];

return;


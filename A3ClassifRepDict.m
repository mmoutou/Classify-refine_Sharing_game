function A3 = A3ClassifRepDict(allP, HIVal, SIVal)
% A3ClassifRepDict - A map for probability of actual fair return.
% midPFair is optional argument, to avoid re-calculating it if already available.
%
%   allP must include, for use here: wH, wS, w0, resNRepo, retLevN .

wH = allP.wH;               wS = allP.wS;                w0 = allP.w0; 
resNRepo = allP.resNRepo;   resNHub = allP.resNHub;
retLevN = allP.retLevN;

Ioff = 0.5;
try 
    if isempty(retLevN); retLevN=3; end
catch
    retLevN = 3;
end


% A{3} : Actual fair ret. does NOT depend only on trueHI and trueSI, in that
%        the initial and final states give indifferent returns
A3 = zeros(retLevN,    resNHub,  resNHub, resNRepo+2); 
A3(3,:,:,resNRepo+1:resNRepo+2) = 1;  % (3,... is indifferent return
for kS = 1:resNHub    
  SI = SIVal(kS);  
  for kH = 1:resNHub
        HI = HIVal(kH);
        % Py of unfair return
        A3(1,kH,kS,1:resNRepo) = 1- 1/(1+exp(w0+(HI-Ioff)*wH + (SI-Ioff)*wS));         
  end
end
A3(2,:,:,:) = 1-(A3(1,:,:,:)+A3(3,:,:,:));  % Py of fair return
% NB we don't need separate noisification for A3

return;


function AHub = AHubClassSerDict(allP, HIVal, SIVal)
% AHubClassRepDict - A map for probability of actual fair return.
%
%   allP must include, for use here: wH, wS, w0, resNHub, retLevN .
%   HIVal and SIVal best reserved for generative process, otherwise omit and use 
%   defaults as per allP.HIv0 , allP.SIv0

wH = allP.wH;               wS = allP.wS;                w0 = allP.w0; 
resNHub = allP.resNHub;
retLevN = allP.retLevN-1;

Ioff = 0.5;
try 
    if isempty(retLevN); retLevN=2; end
catch
    retLevN = 2;
end
try
    SIVal;
catch
    SIVal = allP.SIv0; % [0.25 0.75];
end
try
    HIVal;
catch
    HIVal = allP.HIv0; % [0.25 0.75];
end
% AHub : Actual fair ret. depends only on trueHI and trueSI
AHub = zeros(retLevN,  resNHub,  resNHub); 

for kS = 1:resNHub    
  SI = SIVal(kS);  
  for kH = 1:resNHub
        HI = HIVal(kH);
        % Py of unfair return
        AHub(1,kH,kS) = 1- 1/(1+exp(w0+(HI-Ioff)*wH + (SI-Ioff)*wS));         
  end
end
AHub(2,:,:) = 1-(AHub(1,:,:)); % +AHub(3,:,:));  % Py of fair return

return;


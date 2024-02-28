function AHub = AHubSerDict(allP)
% AHubRepDict - A map for probability of actual fair return.
%               Simplest case ...
%
%   allP must include, for use here: wH, wS, w0, resNHub, retLevN .

wH = allP.wH;               wS = allP.wS;                w0 = allP.w0; 
resNHub = allP.resNHub;
retLevN = allP.retLevN-1;
midGridHub = allP.midGridHub;

Ioff = 0.5;

% AHub : Actual fair ret. depends only on trueHI and trueSI
AHub = zeros(retLevN,  resNHub,  resNHub); 

for kS = 1:resNHub  
  SI = midGridHub(kS);
  for kH = 1:resNHub
        HI = midGridHub(kH);
        % Py of unfair return
        AHub(1,kH,kS) = 1- 1/(1+exp(w0+(HI-Ioff)*wH + (SI-Ioff)*wS));         
  end
end
AHub(2,:,:) = 1-(AHub(1,:,:)); % +AHub(3,:,:));  % Py of fair return

return;


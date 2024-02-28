function [mdpHI, mdpSI] = hub2spoke_DAttr_vii(mdpHI,mdpSI, mdpHubElement, allP)
% function [mdpHI, mdpSI] = hub2spoke_dAttr(mdpHI,mdpSI, mdpHubElement, allP, pFHS, attrV, iLo, iHi, iHS)
%          update D in the attribution MDPs according to the current mdpHub
%          for spm_mdp_L_vii . 
% For this, we need to know the posterior d and a for the hub MDP, to see what attributes the new
% beliefs about a and d amount to. To derive what the 'hi' and 'lo' states
% now amount to, we'll  use the portion of aHub that maps 
% from hi-lo HI and SI states to returns.
  
  resNRepo = allP.resNRepo;
  midGridHub = allP.midGridHub;
  noiseFloor = allP.noiseFloor;
  
  % Prior beliefs about initial states: d map ------------------------------------ 
  dHubHI =  mdpHubElement.d{1}; % Note these are after learning d, which is ...
  dHubSI =  mdpHubElement.d{2}; % ... what we need for the input of the late part of the trial here.
  
  % total evidence and means:
  evHI = sum(dHubHI);            pHubHI = dHubHI / evHI;
  evSI = sum(dHubSI);            pHubSI = dHubSI / evSI;
  mHI = sum(pHubHI .* midGridHub);   
  mSI = sum(pHubSI .* midGridHub);
  
  % Approximation for HI attribution MDP: 
  alphaHI = mHI * evHI; 
  betaHI = evHI - alphaHI;  
  cbd = [0 betacdf( (1:(resNRepo-1))/resNRepo, alphaHI, betaHI ) 1];
  pmf = cbd(2:end) - cbd(1:resNRepo) + noiseFloor; 
  mdpHI.D{1} = (pmf' / sum(pmf)) * sum(dHubHI);
  
  % Similar for SI attribution MDP: 
  alphaSI = mSI * evSI; 
  betaSI = evSI - alphaSI;  
  cbd = [0 betacdf( (1:(resNRepo-1))/resNRepo, alphaSI, betaSI ) 1];
  pmf = cbd(2:end) - cbd(1:resNRepo) + noiseFloor; 
  mdpSI.D{1} = (pmf' / sum(pmf)) * sum(dHubSI);
 
  % F2 are 'report of level of Intent', certain to be at their initial state:
  mdpHI.D{2} = [ones(resNRepo,1); 512];
  mdpSI.D{2} = [ones(resNRepo,1); 512];
  
return


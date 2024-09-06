function [mdpHI, mdpSI] = hub2spoke_dAttr(mdpHI,mdpSI, mdpHubElement, allP)
% function [mdpHI, mdpSI] = hub2spoke_dAttr(mdpHI,mdpSI, mdpHubElement, allP, pFHS, attrV, iLo, iHi, iHS)
%          update d in the attribution MDPs according to the current mdpHub
% For this, we need to know the posterior d and a for the hub MDP, to see what attributes the new
% beliefs about a and d amount to. To derive what the 'hi' and 'lo' states
% now amount to, we'll  use the portion of aHub that maps 
% from hi-lo HI and SI states to returns.
  
  resNRepo = allP.resNRepo;
  retLevN = allP.retLevN; 
  noiseFloor = allP.noiseFloor;
  pFHS = allP.pFHiLoHS.pFHS;
  attrV = allP.pFHiLoHS.attrV;
  iLo = allP.pFHiLoHS.iLo;                    iLen = length(iLo)^2; 
  iHi = allP.pFHiLoHS.iHi;
  iHS = allP.pFHiLoHS.iHS; 

  aHISI2Fair = mdpHubElement.a{1}(1:(retLevN-1),:,:);    
  eLoHi = sum(aHISI2Fair,1); eLoHi = eLoHi(:);  % Total evidence for each beta distro pUnfair--pFair
  mLoHi = aHISI2Fair(2,:,:); mLoHi = mLoHi(:) ./ eLoHi;   % means. NB a{3} row 2 is fair / good.
  % Square-deviates of means from the corresponding probabilities produced by possible HI, SI,
  % weighed by the amount of evidence for each such mean:
  wDist2 = (pFHS - repmat(mLoHi',iLen,1)).^2 ;          % Square deviates 
  wDist2 = sum(wDist2 .* repmat(eLoHi',iLen,1) , 2);    % weigh by evidence and sum
  [~, iWDmin] = min(wDist2);  % find the most representative set HILo-SILo, HIHi-SILo, HILo,SIHi, HIHi,SIHi  
  % ########  Crucial 'correct' Likert indices for Hi and Lo states ########### 
  iHLo = iLo( iHS( iWDmin, 1 ) );        iHHi = iHi( iHS( iWDmin, 1 ) );
  iSLo = iLo( iHS( iWDmin, 2 ) );        iSHi = iHi( iHS( iWDmin, 2 ) );
  % The actual values of Hlo,Hhi etc:
  mdpHI.attrVLoHi = [attrV(iHLo)  attrV(iHHi)]';
  mdpSI.attrVLoHi = [attrV(iSLo)  attrV(iSHi)]';
  
  % Prior beliefs about initial states: d map ------------------------------------ 
  dHubHI =  mdpHubElement.d{1}; % Note these are after learning d, which is ...
  dHubSI =  mdpHubElement.d{2}; % ... what we need for the input of the late part of the trial here.
  % Approximate effective counts, treating Hlo etc. as if proportions. 
  % This is handwavey, I hope it will do ... :
  
  % For HI attribution MDP: 
  alphaHI = sum( mdpHI.attrVLoHi .* dHubHI ); 
  betaHI = sum(dHubHI) - alphaHI;  
  cbd = [0 betacdf( (1:(resNRepo-1))/resNRepo, alphaHI, betaHI ) 1];
  pmf = cbd(2:end) - cbd(1:resNRepo) + noiseFloor; 
  mdpHI.d{1} = (pmf' / sum(pmf)) * sum(dHubHI);
  
  % Similar for SI attribution MDP: 
  alphaSI = sum( mdpSI.attrVLoHi .* dHubSI ); 
  betaSI = sum(dHubSI) - alphaSI;  
  cbd = [0 betacdf( (1:(resNRepo-1))/resNRepo, alphaSI, betaSI ) 1];
  pmf = cbd(2:end) - cbd(1:resNRepo) + noiseFloor; 
  mdpSI.d{1} = (pmf' / sum(pmf)) * sum(dHubSI);
 
  % F2 are 'report of level of Intent', certain to be at their initial state:
  mdpHI.d{2} = [ones(resNRepo,1); 512];
  mdpSI.d{2} = [ones(resNRepo,1); 512];
  
return


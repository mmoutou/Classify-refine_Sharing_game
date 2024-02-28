function fairRep_D = hub2spoke_DFair_vii(hub_d, allP )
% function mdpFairRep = hub2spoke_DFair_vii(mdpFairRep,hub_a1, allP )
% Provides normalised prob. D : fairRep_D{1} = upmf' / sum(upmf);  
%          Estimate the D map for reporting expected fairness, based on the
% effective counts in the a map of the hub. That is, we assume that the cells of
% the a map summarise all effective observations that the pt considers for
% this particular interaction.
%      hub_a1 is the appropriate trial mdpHub(trial).a{1}

  resNRepo = allP.resNRepo;
  resNHub = allP.resNHub;
  noiseFloor = allP.noiseFloor;
  Ioff = allP.Ioff;
  w0 = allP.w0;
  wS = allP.wS;
  wH = allP.wH;
  midGridHub = allP.midGridHub ;

  % the uncertainty in the state does not contribute to the
  % final uncertainty !!! :
  ph = hub_d{1}/ sum(hub_d{1}) ;
  ps = hub_d{2}/ sum(hub_d{2}) ;
  pHS = ph * ps';
  
  pFhs = zeros(resNHub,resNHub);
  for kS = 1:resNHub
      SI = midGridHub(kS);
      for kH = 1:resNHub
          HI = midGridHub(kH);
          pFhs(kH,kS)= 1/(1+exp(w0+(HI-Ioff)*wH + (SI-Ioff)*wS));
      end
  end 
  
  V = pFhs .* pHS;
  mV = sum(V(:));        % expected fair return
  
  % An approximation to the effective number of observations that have
  % contributed: If there really was N observations, each cell would have
  %  ph_i ps_j N = ph_i ps_j sqrt(Ns*Nh) , as we'd have N=Ns=Nj . 
  % So, consistent quantity would be 
  nEff = sqrt( sum(hub_d{1}) * sum(hub_d{2}) ) ;
  
  % Moment-matched beta distro with the mean and evidenc of the 
  % mean estimate: 
  A = mV * nEff; % was: [A, B] = betaMS2ab (M, SE);
  B = nEff - A;
  cbd = [0 betacdf( (1:(resNRepo-1))/resNRepo, A, B) 1];
  upmf = (cbd(2:end) - cbd(1:resNRepo)) + noiseFloor ;     % add a bit of noise
  fairRep_D{1} = upmf' / sum(upmf);                    % normalise

  % F2 are 'report of level of Intent', certain to be at their initial state:
  fairRep_D{2} = [ones(resNRepo,1); 512]/(resNRepo+512);

return


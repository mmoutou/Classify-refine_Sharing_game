function AOne = AClassifAttrRepDict(allP,kLoHi)
% AClassifAttrRepDict - 
This is a draft attempting to map directly from hub to attribution mdps -
works but produces polymodal attribution mdfs which isn't what we want.

How well one attribution (e.g. HI) is to be reported - written
%    for separate-attribution-MDP spokes, and for version where the reporting of 
%    attributions may happen on a different granularity than the classification of
%    underlying states, cf. serialDictator08
% 
% allP must include, for use here: desBias, Ucor, corLevN,
%                    resNHub (number of e.g. HI states), resNRepo (resolution of
%                    report states)

desBias = allP.desBias; Ucor=allP.Ucor; corLevN=allP.corLevN;
resNHub = allP.resNHub; 
resNRepo = allP.resNRepo;
% pSuccMin = 0.05;           pSuccMax = 0.95;           % See below
if corLevN == 3
      corGrid = [0.1, 0.6, 0.9];
else
      corGrid = 1/(2*corLevN)+(0:(corLevN-1))/corLevN;
end

tiny = 1e-16;    % probabilities of this order are zero for all intends and purposes.
%% Arguments that may be provided, ommitted altogether, or defaulted by [].

% If a particular Likert-scale-like report state denotes a bin centered at attribution fraction repFr, 
% while the hub state has true attribution hubFr, the param of the noisyBino will be
% pSucc = pSuccMax - (pSuccMax - pSuccMin)* abs(repFr - hubFr); 


%% Calculate AOne : How well Intent is to be reported.
%             (row)         (col)         (page)                  
%         report quality   trueIntent   IntentReport   
AOne =    zeros(corLevN+1,   resNHub,    resNRepo+1);   % how well HI or SI reported
% At initial report state, all other state components 
% return indifferent/neutral, i.e. probability 1 for corLevN+1 :
AOne(corLevN+1,:,resNRepo+1  ) = 1;    

% Now calc. the resNHub cols for time 2 on the basis of the trueIntent, 
% e.g. page A3(3, HIlev, SILev, 4, SIrep ) would be high if Intent-lev-to-report=4,
% *shifted by the social desirability bias* ... 
%    To have several levels to account for big or small misattributions. Say we
% call diff of 0 OK, up to 1 lev wrong OKish, the rest just wrong.
for kIntRep = 1:resNRepo         %  Reported Intent level index
     for kIntTru = 1:resNHub     % HI/SItrue (as believed!) index, from lowest to highest value
                                 % of the attribution in question, e.g. 1 for low-HI 2 for high-HI.
                                 % or 1 for low-SI, 2 for medium-SI, 3 for high-SI, or such.
              % first, the unbiased distibution over correctness levels:
              err = abs(kIntRep - kLoHi(kIntTru));  % error, in terms of levels
              if err == 0
                  corPar = corGrid(corLevN);     % correct case - 'corPar' is prob. of 'correct' outcome
              elseif err == 1
                  corPar = corGrid(corLevN-1);   % approx. OK
              else
                  corPar = corGrid(1);           % 'wrong'
              end
              % Now shift it acc. to desirability bias. If desBias is 0, leave corPar alone.
              if desBias > tiny         % +ve desir bias is about under-reporting HI, SI
                  if kIntRep >= kLoHi(kIntTru)  % IntRep is high enough for pmf to be shifted 
                                                % towards wrong
                     corPar = (1-desBias)*corPar;
                  else                % IntRep is on the low side, so shift towards correct
                     corPar = corPar+desBias*(1-corPar); 
                  end
              elseif desBias < -tiny  % -ve desir bias is about over-reporting HI, SI.
                  if kIntRep > kLoHi(kIntTru)   % IntRep is high enough for pmf to be shifted 
                                                % towards correct
                     corPar =  corPar-desBias*(1-corPar);
                  else                % IntRep is on the low side, so shift towards wrong
                     corPar = (1+desBias)*corPar; 
                  end
              end  % end if statement that shifts corPar according to desirability bias desp.
              
              % Load pmf over correctness levels into AOne 
              %
              AOne(1:corLevN,kIntTru,kIntRep) = noisyBino( corPar, Ucor, corLevN ); 
              %
     end  % end loop over all the kIntTrue for this kIntRep
end  % end loop over all kIntRep

return;


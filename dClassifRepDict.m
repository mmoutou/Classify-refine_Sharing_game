function d1 = dClassifRepDict(allP)
%DLEARNREPDICT form d1 map - priors over HI, SI and reporting fair return state 
%   ... on the  basis of central tendency of belief about Harmful intent, pH0,
%   pS0, and the amount of 'notional evidence' for these, dInit.
%   This fn. is for the repeated dictator gain
%
%   allP must include, for use here: pH0, pS0, dInit, resNHub, resNRepo. May include
%        Upop, the dispersion param. for the population distro.

pH0 = allP.pH0; pS0 = allP.pS0;   
try Upop = allP.Upop;  catch Upop = 1;  end
dInitH = allP.dInit;
try   % early versions of the code did not have dInitRat , the dInit S/H ratio.
   dInitS = allP.dInitRat * dInitH; 
catch
   dInitS = dInitH;
end
resNHub=allP.resNHub;   

% Factor 1 is priors over 'true HI' state 
p = (noisyBino(pH0,Upop,resNHub)+0.05)';       d1{1} = dInitH*p/sum(p) ; 
% Factor 2 is priors over 'true SI' state
p = (noisyBino(pS0,Upop,resNHub)+0.05)';       d1{2} = dInitS*p/sum(p) ; 
% F3 is 'report of level of fair return', certain to be at its initial state:
d1{3} = [ones(allP.resNRepo,1); 512; 1];


return;


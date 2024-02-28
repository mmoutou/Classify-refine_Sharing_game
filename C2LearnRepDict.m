function C2 = C2LearnRepDict(allP)
% C2LearnRepDict - factor 2 is: Report of chance of fair return is correct or not.
% desCorr: currency of desirability to be correct
%
%   allP must include, for use here: desCorr, corLevN .

desCorr = allP.desCorr;   corLevN = allP.corLevN; Tsteps1 = allP.Tsteps1; 

weights =  [-3,-1,4,0]'; % 4th row is indifferent reference
         
% Outcome factor 2 is: Report of chance of fair return is correct or not.
% Only at the second timestep will it be important
% important to make the right report, and NOT to make no report:
C2  = zeros(corLevN+1,Tsteps1);        % most are neither good nor bad, but ...
% ... timepoint 2 is non-trivial :
C2(:,2) = desCorr*weights; 

return;


function COne = COneAttrRepDict(allP)
% COneAttrRepDict - factor 3 of Attributing MDP is desirability to report what one wants to.
%   This fn. applies equally to all intent factors in their corresponding mini-mdps.
% desCorr: currency of desirability to be correct
%
%   allP must include, for use here: desCorr, corLevN, Tsteps2

desCorr = allP.desCorr;  corLevN = allP.corLevN;  Tsteps2 = allP.Tsteps2;
weights = [-4,-2,3,-1]'; % 4th row is indifferent reference

% Only at the second timestep will it be important
%  to make the right report, and NOT to make no report:
COne  = zeros(corLevN+1,Tsteps2);   % most are neither good nor bad, but ...
% ... timepoint 2 is non-trivial :
COne(:,2)  = desCorr*weights; % 4th row is indifferent reference

return;


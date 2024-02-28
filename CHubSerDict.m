function C = CHubSerDict(allP)
% CHubRepDict - gives: desirability of fair etc. split
% desFair: currency of desirability for fair split
%
%   allP must include, for use here: desFair

weights = [-1,1]'  ; % rows unfair, fair (only given at one timestep)

C  =  [ [0 0]' weights]*allP.desFair; 

return;


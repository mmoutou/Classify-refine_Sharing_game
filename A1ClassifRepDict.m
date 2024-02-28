function A1 = A1ClassifRepDict(allP)
% A1ClassifRepDict - simply 'where we are' Informational, Report stage
% for the version of repeated dictator where the hub is about classification into
% lo / (med / ) hi rather than many levels, but then there's learning

resNRepo = allP.resNRepo;     % Resolution of reporting HI / SI, e.g. 6
resNHub  = allP.resNHub;      % Resol. of classifying pls HI / SI, e.g. lo - mid - hi only

% A1 : Informational, Report stage
%             (row)              (col)    (page)   (issue)           
%             outcome-to-report  trueHI   trueSI   [p(fair)_level,init,final]  
A1 = zeros(   resNRepo+1,        resNHub, resNHub, resNRepo+2); 
% Report outcomes from initial report state, resolN+1, as well as the
% finall report state, resolN+2, are all 'null', aka noReportMade, i.e. resolN+1:
A1(resNRepo+1,:,:,resNRepo+1:resNRepo+2) = 1;
% Report outcomes from valid report states are as themselves:
for toRep = 1:resNRepo
    A1(toRep,:,:,toRep) = 1;
end

return;


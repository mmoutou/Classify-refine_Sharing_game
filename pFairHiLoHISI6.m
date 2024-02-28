function pFHiLoHS  = pFairHiLoHISI6(wVec,resNRepo,Ioff,attribVec)
% pFairHiLoHISI6  return pFair for HILo-SILo, HIHi-SILo, HILo-SIHi, HIHi-SIHi
%                 combinations of binary attributes Lo or Hi, of 2 attributes, 
%                 HI and SI.
% pFHiLoHS provides:
%   pFHS         gives the four-probability vectors
%   attrV        gives the 6 values of the attribute that each HI and SI state may be
%   iLo and iHi  give the order of the indices (15 each, so pFHS is 15^2 combination x 4 pFairs)
%
% test with: pFHiLoHS = pFairHiLoHISI6([0.5,4,6])
if  resNRepo ~= 6
    error('resNRepo not 6, so pFairHiLoHISI6 not appropriate');
end
try
    attrV = attribVec;  
catch
    % Let's only consider the following values of Attributions, as generated below:
    % [0.0833    0.2500    0.4167    0.5833    0.7500    0.9167 
    attrV = 1/(2*resNRepo)+(0:(resNRepo-1))/resNRepo;  % grid of midpoints at specified resolution
end   
if length(attrV) ~= 6; error('length(attibVec) must be 6'); end
try
    Ioff;
catch
    Ioff = 0.5;
end
w0 = wVec(1); wH = wVec(2); wS = wVec(3);   % params for logistic prob fn below
likLen = 4;   % length of [pHILoSILo, pHIHiSILo, pHILoSIHi, pHIHiSIHi]

%% Formulate possible pairs of lo - hi attribution values
% Now each possible lo-hi pairs with distinct endpoints will
% have indices over these as follows, so e.g. the 2nd iLo and iHi will
% define attribution states  0.0833,  0.4167 
iLo = [ 1 1 1 1 1   2 2 2 2  3 3 3   4 4   5];
iHi = [ 2 3 4 5 6   3 4 5 6  4 5 6   5 6   6];
iLen = length(iLo); 

%% Calculate the p(fair) corresponding to each possible combination of attributes,
%  HILo-SILo, HIHi-SILo, HILo,SIHi, HIHi,SIHi

pFHS = zeros(iLen^2,likLen); 
iHS  = zeros(iLen^2,2); 

cnt = 0;
for iS = 1:iLen
    for iH = 1:iLen
        cnt = cnt+1;
        iHS(cnt,:)= [iH,iS];        % Convenient copy of indices for H and S
                                    % within iLo, iHi.
        HILo = attrV( iLo( iH ) );
        HIHi = attrV( iHi( iH ) );
        SILo = attrV( iLo( iS ) );
        SIHi = attrV( iHi( iS ) );
        
        %  store in order: HILo-SILo, HIHi-SILo, HILo-SIHi, HIHi-SIHi
        pFHS(cnt,1) = 1/(1+exp(w0+(HILo-Ioff)*wH + (SILo-Ioff)*wS));
        pFHS(cnt,2) = 1/(1+exp(w0+(HIHi-Ioff)*wH + (SILo-Ioff)*wS));
        pFHS(cnt,3) = 1/(1+exp(w0+(HILo-Ioff)*wH + (SIHi-Ioff)*wS));
        pFHS(cnt,4) = 1/(1+exp(w0+(HIHi-Ioff)*wH + (SIHi-Ioff)*wS));
    end
end

pFHiLoHS.pFHS = pFHS; 
pFHiLoHS.attrV = attrV;
pFHiLoHS.iLo = iLo;
pFHiLoHS.iHi = iHi;
pFHiLoHS.iHS = iHS; 

return;

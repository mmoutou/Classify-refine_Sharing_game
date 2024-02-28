function pFHiLoHS  = pFairHiLoAttr(wVec,grN,Ioff,attribVec)
% pFairHiLoAttr return pFair for HILo-SILo, HIHi-SILo, HILo-SIHi, HIHi-SIHi
%                 combinations of binary attributes Lo or Hi, of 2 attributes, 
%                 HI and SI, but in the ...Attr version the HI and SI intevals have
%                 the same pair of hi - lo values.
% pFHiLoHS provides:
%   pFHS         gives the four-probability vectors
%   attrV        gives the 6 values of the attribute that each HI and SI state may be
%   iLo and iHi  give the order of the indices (15 each, so pFHS is 15^2 combination x 4 pFairs)
%
% test with: pFHiLoHS = pFairHiLoAttr([0.5,4,6])
try
    grN;
catch
    grN=10;
end
try
    attrV = attribVec; 
    if length(attrV) ~= grN; error('attrV must have grN elements'); end
catch
    % Let's only consider the following values of Attributions, as generated below:
    % [0.0833    0.2500    0.4167    0.5833    0.7500    0.9167 
    attrV = 1/(2*grN)+(0:(grN-1))/grN;  % grid of midpoints at specified resolution
end   
try
    Ioff;
catch
    Ioff = 0.5;
end

w0 = wVec(1); wH = wVec(2); wS = wVec(3);   % params for logistic prob fn below
likLen = 4;   % length of [pHILoSILo, pHIHiSILo, pHILoSIHi, pHIHiSIHi]

%% Formulate possible pairs of lo - hi attribution values
% Now each possible lo-hi pairs with distinct endpoints will
% have indices over these as follows: 
iLen = (grN*(grN-1))/2;   % worked out algebraically ...
iLo = zeros(1,iLen);      iHi = iLo;     
cnt = 1;
for lo = 1:(grN-1)
    for hi = (lo+1):grN
        iLo(cnt) = lo;
        iHi(cnt) = hi;
        cnt = cnt+1;
    end
end

%% Calculate the p(fair) corresponding to each possible combination of attributes,
%  HILo-SILo, HIHi-SILo, HILo,SIHi, HIHi,SIHi

pFHS = zeros(iLen,likLen); 
iHS  = [1:iLen; 1:iLen]';   % Convenient copy of indices for H and S  within iLo, iHi.
                           % (in this version of pFairHiLo ... they are the same.
for cnt = 1:iLen

    HILo = attrV( iLo( cnt ) );
    HIHi = attrV( iHi( cnt ) );
    SILo = HILo;
    SIHi = HIHi;
        
    %  store in order: HILo-SILo, HIHi-SILo, HILo-SIHi, HIHi-SIHi
    pFHS(cnt,1) = 1/(1+exp(w0+(HILo-Ioff)*wH + (SILo-Ioff)*wS));
    pFHS(cnt,2) = 1/(1+exp(w0+(HIHi-Ioff)*wH + (SILo-Ioff)*wS));
    pFHS(cnt,3) = 1/(1+exp(w0+(HILo-Ioff)*wH + (SIHi-Ioff)*wS));
    pFHS(cnt,4) = 1/(1+exp(w0+(HIHi-Ioff)*wH + (SIHi-Ioff)*wS));
end

pFHiLoHS.pFHS = pFHS; 
pFHiLoHS.attrV = attrV;
pFHiLoHS.iLo = iLo;
pFHiLoHS.iHi = iHi;
pFHiLoHS.iHS = iHS; 

return;

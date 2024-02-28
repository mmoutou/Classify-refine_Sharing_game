function pPosHiLoAttr  = pPosHiLoAttr(wVec,grN,Ioff,attribVec)
% pFairHiLoAttr return p positive outcome for AttrLo, AttrHi binary levels of an 
%                 an attribute Attr, 
% pPosHiLoAttr provides:
%   pPosAttr     gives the probability vectors
%   attrV        gives the 6 values of the attribute 
%   iAttr        gives the order of the indices (15 each, so pFHS is 15^2 combination x 4 pFairs)
%
% test with: pPosHiLoAttr  = pPosHiLoAttr(([0.5,4,6])
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

w0 = wVec(1); wAttr = wVec(2);    % params for logistic prob fn below
likLen = 2;   % length of [pAttrLo, pAttrHi]

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
%  AttrLo, AttrHi

pPosAttr = zeros(iLen,likLen); 
iAttr  = [1:iLen; 1:iLen]';   % Convenient copy of indices for Attr  within iLo, iHi.

for cnt = 1:iLen

    AttrLo = attrV( iLo( cnt ) );
    AttrHi = attrV( iHi( cnt ) );
        
    %  store in order: AttrLo, AttrHi
    pPosAttr(cnt,1) = 1/(1+exp(w0+(AttrLo-Ioff)*wAttr ));
    pPosAttr(cnt,2) = 1/(1+exp(w0+(AttrHi-Ioff)*wAttr ));
end

pPosHiLoAttr.pPosAttr = pPosAttr; 
pPosHiLoAttr.attrV = attrV;
pPosHiLoAttr.iLo = iLo;
pPosHiLoAttr.iHi = iHi;
pPosHiLoAttr.iAttr = iAttr; 

return;

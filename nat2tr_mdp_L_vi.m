function trP = nat2tr_mdp_L_vi(P)
% nat2tr_mdp_L_ii - to bring param for mdp_L_ii from native to 
%                  transformed space. natP is a structure w named fields.
%__________________________________________________________________________

% testing = 1;   % for debugging etc.

% Complicated if statement to transform inputted parameters 
field = fieldnames(P);  
for i = 1:length(field)
    % first, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if strcmp(field{i},'initEvRat')
        trP.initEvRat = log(P.initEvRat);   
    elseif strcmp(field{i},'dInitEv')
        trP.dInitEv = log(P.dInitEv);   
    elseif strcmp(field{i},'aInitEv')
        trP.aInitEv = log(P.aInitEv);   
    elseif strcmp(field{i},'alphaPrec')
        trP.alphaPrec = log(P.alphaPrec);        
    elseif strcmp(field{i},'wH')
        trP.wH = log(P.wH);      
    elseif strcmp(field{i},'wS')
        trP.wS = log(P.wS);    
    elseif strcmp(field{i},'Ucor')
        trP.Ucor = log(P.Ucor);   
    % logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'pH0')
        pNat = P.pH0;       trP.pH0 = log(pNat/(1-pNat)); 
    elseif strcmp(field{i},'pS0')
        pNat = P.pS0;       trP.pS0 = log(pNat/(1-pNat));  
    elseif strcmp(field{i},'lrnR')  
        pNat = P.lrnR;       trP.lrnR = log(pNat/(1-pNat));  
    elseif strcmp(field{i},'mem')  
        pNat = P.mem;       trP.mem = log(pNat/(1-pNat));  
    % and scaled-logit transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'desBias')
        pNat =  P.desBias;  trP.desBias = log((1+pNat)/(1-pNat));
    % Assume every thing else Untransformed ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
       trP.(field{i}) = P.(field{i});  
   end
end
  

return;

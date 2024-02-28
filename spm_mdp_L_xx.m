function fitM = spm_mdp_L_xx( trP, ModStruc, Inp, Resp, details)
% 
%  trP is the a parameters structure with all the fields, within and
%  between block, needed. Or alternatively a vector of the same, in order.
%  Inp and Resp is all data for one pt from e.t. expD09{ptN}
%  ModStruc is the MDP model structure
%  Derived from, and using, spm_mdp_L_vi- BUT DOES NOT CATER FOR PRIOR OVER PARAMS.


% BASIC PIPELINE:
%     . serialDictator* -> formulate and practice MPD based model
%                  and assemble MAP structure / DCM structure for model-fitting to be 
%                  used within EstimParAcIn below (in this fn).
% --> . spm_mdp_L_*    ----> Likelihood fn for all the parameters to be studied
%     . spm_mdp_L_*_wrap* -> Likelihood function wrapper for subsets or other variants of above
%     . [experiment]Dat4AcIn*   -> bring exprerimental data to format of generative model
%            e.g. ...\Dropbox\BASOR\AcInSOR\AcInRepeatedDictator\attrssrib\attrSSRIbDat4AcIn09.m
%     [ . [expt]Grid[number][let]_[roman]_[pts block}  -> Grid fit INCL learning from one block to the next
%     { . EstimParAcIn*   -> Use fmincon (or spm_nlsi_Newton...) to fit model structure
%     { . [experiment]Fit[number][letter]_[latin]_[number], e.g. attrssriaFit07a_iv_01 to fit 1st batch of data,
%               with param combination 'a', using likelihood function iv, gen model 07.
%     . mergeExpFits[number][letter]_[number] --> produce nice csv out of multiple fits that will have been 
%               done in parallel, including descriptive data.
%  then, to avoid local minima, if a number of models has been run:
% 

try 
    details;     % details > 0 enriches mL, >= 2 plots stuff.
catch
    details=0;                     % default if actions given is to just output the log lik
end   

% Headings for params, in exact order specified in the model structure 
% below, i.e. modStruc.indexP - then block-level params, then fitting measures.
%          1     2      3          4           5        6          7    8     9    10     11    12    
parHd = {'pH0','pS0','dInitEv','aInitEv','initEvRat','alphaPrec','mem','wH' ,'wS','w0','othLR','pocB'};
totBlN = length(Inp); 
% 'reminder' of 'Likert scale' type bins for reported beliefs of the expected pH0 and pS0:
resNRepo = ModStruc.allP.resNRepo; 
midGrid = 1/(2*resNRepo)+(0:(resNRepo-1))/resNRepo;

%% Restore parameter argument as structure
% REM trP AND indexP have ONLY the 'index' params, e.g. to be fitted! 
%              Not the ones needed for additional specification / fixed ones! 
%              To see all the params, incld he fixed ones, list modStruc.allP .
% First, bring P from its inputted form, such as a vector, to the form that 
% the param structure modStruc.indexP also has :
if ~isstruct(trP); trP = spm_unvec(trP,ModStruc.indexP); end
ModStruc.indexP = trP;
if ~isempty(ModStruc.priPar)
    error('Do not attempt to pass priors to spm_mdp_L_xx ! ModStruc.priPar was not empty.');
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Aggregate over all blocks (computer partners) 
%  Aggregate posterior density for this pt and parameter combination    ~~~~~~~~~
%                              over blocks.                              ~~~~~~~~~
% First, working copies that can be shifted by other-to-other learning and biases
% during the loop over partners:
pH0Tr = trP.pH0;               pS0Tr = trP.pS0;
pH0   = 1/(1+exp(-pH0Tr));     pS0   = 1/(1+exp(-pS0Tr)); 
dInitEv = exp(trP.dInitEv); 
aInitEv = exp(trP.aInitEv); 
% And for convenience (these will not be shifted):
lambda = 1/(1+exp(-trP.othLR));   % native space LR
pocB  = trP.pocB;

% Iterate over all blocks and accumulate log lik etc.
for blockN = 1:totBlN

	%%  Specific Inputs and Responses to each experimental block:
	ModStruc.Inp =  Inp(blockN);
	ModStruc.Resp = Resp(blockN);

                     
    %% A. Modify according to learning through blocks
    if blockN > 1     % i.e. there's been feedback to learn from

        % Apply learning to native, not transformed, space (approximation checked algebraically, to a degree :) )
        % debug line: if blockN ==2; warning('pH0 update:'); disp([pH0 pH0 + lambda*DH0]);  end
        pH0 = pH0 + lambda*DH0; 
        if pH0 > 0.999; pH0 = 0.999; elseif pH0 < 0.001; pH0 = 0.001; end  % Correct 'road exits'!
        pS0 = pS0 + lambda*DS0;
        if pS0 > 0.999; pS0 = 0.999; elseif pS0 < 0.001; pS0 = 0.001; end  % Correct 'road exits'!
                    
        % Will be recorded in transformed space in the structure to be passed to 
        % likelihood fn. after applying bias, below.

        % notional evidence increased as per d <- d(1-lrnRate)+ 1 .
        % We don't bother adjusting d for H and S separately. NB dInitEv (not a member of
        % structure) has been transf. to native :
        dInitEv = dInitEv + lambda;  % Not sure of this ... should be OK!
        aInitEv = aInitEv + lambda;  
        if dInitEv < 0.001; dInitEv = 0.001; warning('dInitEv set to 0.001 :( '); end
        if aInitEv < 0.001; aInitEv = 0.001; warning('aInitEv set to 0.001 :( '); end
        ModStruc.indexP.dInitEv = log(dInitEv); 
        ModStruc.indexP.aInitEv = log(aInitEv);                     
    end
    
    %% B. Modify acc. to coallitional bias depending on ethnicity
    if strcmp('whitish',ModStruc.ethn{blockN})
       coalB = 1/pocB; 
    elseif strcmp('poc',ModStruc.ethn{blockN})
       coalB = pocB;
    else 
       error('modInfo.ethn invalid');
    end
    %  ...will give -ve val for whitish, +ve for POC.
    % Bias applied so that higher bias increases pS0 and pH0 for POCs
    % modInfo.mdpStruc.indexP is used for the actual likelihood / post. density 
    % calculation, and the items below from within it are already log-transf.
    % so that e.g. pocB = 2 gives pH0 -> pH0^(1/2), i.e. increased.
    % NB as baseline, only apply coalitional threat to HI bias. Could model-compare w. SI too ... 
        
    % Record in transformed space:
    ModStruc.indexP.pH0 = log(pH0^(1/coalB) /( 1-pH0^(1/coalB) )); 
    ModStruc.indexP.pS0 = log(pS0/(1-pS0)); 
                
    ModStruc.mdpStruc.priPar = [];  % as priors just redefined above.
                  
    % Fit measure for each block:   ====================================================
    try iP = rmfield(ModStruc.indexP,'othLR'); catch iP = ModStruc.indexP; end
    try iP = rmfield(iP,'pocB'); catch ; end

    blFit = spm_mdp_L_vi(iP, ...
                         ModStruc,...
                         Inp(blockN),...
                         Resp(blockN),1);  % last entry is 1 to give
    %   details, but if 0, it gives the MINUS sum LL.
    %   ================================================================================

    %%  Now calculate prediction errors that will be used to update
    % pH0, pS0 and d for the next block. NB these are notated D in the
    % 'reporting spokes' but derive from updated beliefs d of the 'hub'
    % for HI:
    post = blFit.MDPHI(end).D{1}'; post = post/sum(post); 
    DH0 = sum(post .* midGrid) - pH0;
    % for SI:d
    post = blFit.MDPSI(end).D{1}'; post = post/sum(post); 
    DS0 = sum(post .* midGrid) - pS0; 
    % d doesn't need this, in this simplest of formulations
    
    % Store in results stucture(s) for output
    if details
        if blockN == 1   
            fitM.sLL = blFit.sum;
            % Make space for details:
            if totBlN ~= 4
                fitM.blFit = [];
                for k=1:totBlN; fitM.blFit(end+1) = blFit; end
            else
                fitM.blFit = [blFit blFit blFit blFit];
            end
        else
            fitM.sLL = fitM.sLL + blFit.sum;
            fitM.blFit(blockN) = blFit;
        end
    else  % the more common case - just accumulate logLik
        if blockN == 1
            fitM = blFit.sum;
        else
            fitM = fitM + blFit.sum;
        end
    end
  
end   % end loop over blocks

return;  % whole function.





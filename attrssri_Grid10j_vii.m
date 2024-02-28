function optPar = attrssri_Grid10j_vii( block, wave)
%  optPar = attrssribGrid10j_vii( block, wave) :  MAP fitting, with weakly informative priors given by 
%                        scaled beta distributions. ssria means the first wave, and 10j means
%                        that it's the 10 version of the generative function, j for a large set
%                        of parameters to be fitted/considered. 
%                        Fitting being done with simplest learning from one block to the other,
%                        adding admixture of the posterior to the prior.
%                        Estimate over a grid. vii is the version of the likelihood fn used.
% Modify this at will and just store the version used for a bunch of results in the respective dir.

% clear variables;       block = 1;
debugging  = 1;  % affects grid size, block size etc.
grN = 12;       grSc  = 5;        % will try grids of: par0 + (-grN:grN)*(abs(par0)+1)/grSc ; 
if debugging
    warning('running with debug settings'); 
    grN = 1
    grSc = 2
end 
%       1     2      3          4           5        6     7     8    9     10     11    12     13   14    15    16
hd = {'pH0','pS0','dInitEv','initEvRat','alphaPrec','mem','wH' ,'wS','w0','othLR','pocB','LP','LL','AIC','BIC'};

% NB each pt has 4 partners !! So blocksize = 3 means 12 fits, 
blocksize = 10;  if debugging; blocksize = 3; end
todo = (block-1)*blocksize+1 : block*blocksize;
% adjust last block which mops up a bit more:
if strcmp(wave, 'a')
   if (73-todo(end)) < blocksize; todo = [todo,(todo(end)+1):73]; end
elseif strcmp(wave, 'b')
   if (66-todo(end)) < blocksize; todo = [todo,(todo(end)+1):66]; end
else
   error(['wave ' wave '??? ']);
end

cwd = pwd;
wrkDri = ['C:\Users\mmpsy\Dropbox\BASOR\AcInSOR\AcInRepeatedDictator\attrssri' wave '\'];
try
   cd(wrkDri);
catch
   wrkDri = ['/home/mmoutou/Dropbox/BASOR/AcInSOR/AcInRepeatedDictator/attrssri' wave '/'];
   cd(wrkDri);
end

% Loads 'modStruc', which has the key elements of the model structure,
% and expD09 which has the real data, and Pall_template which is a Pall that 
% will be used as a template or scaffold for all the fixed params.
AIDatFile  = ['attrssri' wave 'AcIn09.mat'];            
load(AIDatFile);        % has expD09. Data must be (re-formated) to fit w model Inp, Resp etc.
expDatFile = ['attrssri' wave 'Dat.mat'];               
load(expDatFile);       % This doesn't depend on the generative model.
modStrucFile = 'modStruc10a.mat';              load(modStrucFile);     % This should contain modStruc
datNPerPt = 3 * 4 * 12;  % {data = fairnes prediction, HI, SI ratings} x 4 partners x 12 trials per partner

fitFileName = ['exp_' wave '_Grid10j_' num2str(todo(1)) 'to' num2str(todo(end))];

expGrid = {};
indexPLen = length(fields(modStruc.indexP)); 
par2fitN=(indexPLen)+2;      % We would like to have the choice to fit any of the params in modStuc.indexP
                % but not nessarily all of them. e.g. length(fields(modStruc.indexP)) = 10
                % but only fit 1st 6, so (indexPLen)+2 has all in indexP plus
                % the 7th param is the leaning rate from partner to partner, othLR
                % 8th is the poc bias, pocB.
totParN=indexPLen+2; % total number of params in indexP, plus the two explored here.
totOthN=4; 
% 'reminder' of 'Likert scale' type bins for reported beliefs of the expected pH0 and pS0:
resNRepo = modStruc.allP.resNRepo; 
midGrid = 1/(2*resNRepo)+(0:(resNRepo-1))/resNRepo;
optPar = nan(blocksize, totParN+4);   % for optimised parameters AND 4 measures of model fit,
         % MAPd, LL, AIC, BIC incl the non-MDP params and the ones we won't explore here.

cnt = 0;          % counts all blocks
ptCnt = 0;        % counts participans
for ptN = todo
    ptCnt = ptCnt+1; 
    parOrdGrid = [1 2 8 3 5 11 4 6 7 9 10 1 2 8 3 5 12 7 6 10 4 9 1 2 11];  % in order above. Cycle twice!
    if debugging; parOrdGrid = [7 11 8 11]; end
    parCnt = 0;
    for par2grid = parOrdGrid
        parCnt = parCnt + 1;
        expGrid{ptN,par2grid}.parPDens = nan(3,2*grN+1);  % to work on finding MAP etc. Rows of 
                                                          % param value, fit measure, sum LL.
   
        %% Clear the decks and set stuff common to ALL blocks NOTE NOT ADJUSTED W.R.T. COALITIONAL IDENTITY.
        %  These are meant to be v v weakly informative so that should be OK ...
        % Parameters to fit (the rest remain as within synthD{ptN}.Pall below):
        if cnt > 1; clear mapMod; end
        mapMod.mdpStruc = modStruc;
        % So the othLR has flat prior over unit intrval and othB has approx gaussian.
        % See serialDictator09b for example.
        mapMod.blockParPri = [[1.01,1.01,0,1]' [10,10,-46,46]'];

        %% Starting param values for the (iterative convergence) fit. Some derived from 1st block.
        %  Use the average values of HI and SI as starting points for the priors. NOTE
        %  THESE ARE NOT PRIORS IN THE COG MODEL, ONLY INITIAL VALUES FOR FITTING.
        %  Also master copy of MAP model specification.
        
        if parCnt == 1
           natPinit = tr2nat_mdp_L_vi(modStruc.indexP);  % intialize native space params
                                                         % also for version vii of likelihood.
           % Params 1 and 2 - used inside MDPs :
           blockN = 1; 
           try
              natPinit.pH0 =   (expD{ptN}.partner(blockN).descr(6)-1)/4;   % rescaled from 1-5 to 0-1 
              natPinit.pS0 =   (expD{ptN}.partner(blockN).descr(7)-1)/4;   % rescaled from 1-5 to 0-1 
           catch  
              warning('Could not set initial values from data - BEWARE!');
              natPinit.pH0 = 0.2;    natPinit.pS0 = 0.2;  
           end

           % params 3 to 10 - used inside MDPs :
           natPinit.dInitEv = 1;      % 5: try starting from a strongish value!
           % NOT in version 10: natPinit.aInitEv = 1;   % 5: try starting from a strongish value!
           natPinit.initEvRat  =1;    % boring value ...
           natPinit.alphaPrec =0.5; % 0.5: weak value ...
           natPinit.mem  =0.999;    % no forgetting default
           natPinit.wH  = 6;        % wS and wH added for version attrssriaFit07c_iv_**         
           natPinit.wS  = 6;        % 
           natPinit.w0  = 0.2;      % No Eye Deer if any good ...
           
           % Params 11 to 102, not used inside MDPs but between them:
           othLR = 0.2;                     % NEW FOR attrsri$07b_iv_* : learning rate about others, i.e. between blocks.
           othLRTr = log(othLR/(1-othLR));  % To real line.
           mapMod.mdpStruc.othLRTr = othLRTr;
           pocB = 1.001;                    % Slightly biased value
           pocBTr = log(pocB);              % ... to real line.
           mapMod.mdpStruc.pocBTr = pocBTr; 
        
           % mapMod.indexP is used directly by spm_mdp_L_vi to spm_unvec its 
           % parameter vector argument, whereas fitting routine also uses mapMod.iniTrP, the 
           % initial values for the fit in transformed space :
           mapMod.mdpStruc.indexP = nat2tr_mdp_L_vi(natPinit); 
           fields4MDPs = fieldnames(mapMod.mdpStruc.indexP); 
           refMapMod = mapMod;     % master copy
        end
        
        %% Iterate over the grid
        for iGr = -grN:grN      
            mapMod = refMapMod; % REM mapMod is max a posteriori focused model whose elements
                                % will be modified acc. to learning, coalitional bias etc. and
                                % then be directly passed to likelihood (incl. MAP density) function.

            %% Parameter value acc. to grid, but before adjustment for which block we're in
            if par2grid <= par2fitN-2  % i.e. if it's not the learning-between-others learning rate
                                     % or the pocB
                refp = refMapMod.mdpStruc.indexP.(fields4MDPs{par2grid}); 
                thisp =  refp + iGr * (abs(refp) + 1)/grSc ; 
                mapMod.mdpStruc.indexP.(fields4MDPs{par2grid}) = thisp; 
                parName = [fields4MDPs{par2grid} 'Tr'];
            elseif par2grid == par2fitN-1     % i.e. if stepping the others-lr
                refp = othLRTr;
                thisp =  refp + iGr * (abs(refp) + 1)/grSc;                 
                mapMod.mdpStruc.othLRTr = thisp; 
                parName = 'othLRTr';
            else   % I.e. 
                refp = pocBTr;
                thisp =  refp + iGr * (abs(refp) + 1)/grSc ; 
                mapMod.mdpStruc.pocBTr = thisp; 
                parName = 'pocBTr'; 
            end
            if iGr == -grN && parCnt ==1  % initialise optimal param estimate - elements will be replaced 
                % as we optimise. Last two are for the fit measures. 
               optPar(ptCnt,:)= [spm_vec(mapMod.mdpStruc.indexP)' othLRTr pocBTr nan nan nan nan];
            end
        
            %% Aggregate over all blocks (computer partners) 
            %  Aggregate posterior density for this pt and parameter combination    ~~~~~~~~~
            %                              over blocks.                              ~~~~~~~~~
            % First, working copies that will be shifted by other-to-other learning and biases
            % during the loop over partners:
            pH0Tr = mapMod.mdpStruc.indexP.pH0;        pS0Tr = mapMod.mdpStruc.indexP.pS0;
            pH0   = 1/(1+exp(-pH0Tr));                 pS0   = 1/(1+exp(-pS0Tr)); 
            dInitEv = exp(mapMod.mdpStruc.indexP.dInitEv); 
            % NOT in version 10 : aInitEv = exp(mapMod.mdpStruc.indexP.aInitEv); 

            for blockN = 1:length(expD{ptN}.other)
                cnt = cnt+1;
                
                % display progress:
                t = clock; tstr = num2str(t(2:5)); 
                if cnt == 1; tStart = tstr; end
                if blockN == 1; disp(['Now pt.=' num2str(ptN) ' for ' fitFileName ' at ' tstr]); end

                %%  Specific Inputs and Responses to each experimental block:
                mapMod.Inp =  expD09{ptN}.partner{blockN}.Inp;
                mapMod.Resp = expD09{ptN}.partner{blockN}.Resp;
        
                % For convenience:
                mapMod.other = blockN;
                mapMod.file = expD{ptN}.file;
                mapMod.age = expD09{ptN}.partner{blockN}.age;

                %% modify parameters for this block. NB slighly problematic re. prior densities ...
                %  ... which aren't modified ...
                
                expGrid{ptN,par2grid}.parPDens(1,iGr+grN+1) = thisp;   % make a record
                      
                %% A. Modify according to learning through blocks
                if blockN > 1     % i.e. there's been feedback to learn from
                    lambda = 1/(1+exp(-mapMod.mdpStruc.othLRTr));   % native space LR
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
                    % NOT in mk 10: aInitEv = aInitEv + lambda;  
                    if dInitEv < 0.001; dInitEv = 0.001; warning('dInitEv set to 0.001 :( '); end
                    % NOT in mk 10: if aInitEv < 0.001; aInitEv = 0.001; warning('aInitEv set to 0.001 :( '); end
                    mapMod.mdpStruc.indexP.dInitEv = log(dInitEv); 
                    % NOT in mk 10: mapMod.mdpStruc.indexP.aInitEv = log(aInitEv);                     
                end
                %% B. Modify acc. to coallitional bias depending on ethnicity
                pocB = exp(mapMod.mdpStruc.pocBTr); 
                mapMod.ethn = expD09{ptN}.partner{blockN}.ethn;
                if strcmp('whitish',mapMod.ethn)
                    coalB = 1/pocB; 
                elseif strcmp('poc',mapMod.ethn)
                    coalB = pocB;
                else 
                    error('mapMod.ethn invalid');
                end
                %       will give -ve val for whitish, +ve for POC.
                %  Bias applied so that higher bias increases pS0 and pH0 for POCs
                % mapMod.mdpStruc.indexP is used for the actual likelihood / post. density 
                % calculation, and the items below from within it are already log-transf.
                % so that e.g. pocB = 2 gives pH0 -> pH0^(1/2), i.e. increased.
                pH0 = pH0^(1/coalB); 
                % NB as baseline, only apply coalitional threat to HI bias. Could model-compare w. SI too ... 
                % pS0 = pH0^(1/coalB); 
        
                % Record in transformed space:
                mapMod.mdpStruc.indexP.pH0 = log(pH0/(1-pH0)); 
                mapMod.mdpStruc.indexP.pS0 = log(pS0/(1-pS0)); 
                
                %% ~~~~~~~~~~~ Key Log Lik Density ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                          
                % Copy the MDP building blocks in place:
                mapMod.mdpStruc.iniTrP = mapMod.mdpStruc.indexP;      
                mapMod.mdpStruc.hmmHub = modStruc.hmmHub;
                mapMod.mdpStruc.mdpHI = modStruc.mdpHI; 
                mapMod.mdpStruc.mdpSI = modStruc.mdpSI; 
                mapMod.mdpStruc.allP = modStruc.allP;       
                
                %   ================================================================================
                mL = spm_mdp_L_vii(mapMod.mdpStruc.indexP,mapMod.mdpStruc,...
                     mapMod.Inp,mapMod.Resp,1);  % last entry is 1 to give both
                %   ================================================================================
                   
                % Initiate with sum-log-density over 'external to MDPs' params:
                if blockN ==1
                   % The priors term of the MAP density is only included first time
                   % around. First, the extra params:
                   ptPDens =  sum(dbetasc([othLR,pocB],mapMod.blockParPri(1,:),...
                                        mapMod.blockParPri(2,:),...
                                        mapMod.blockParPri(3,:),mapMod.blockParPri(4,:)));
                   ptSLL = 0; 
                   % Then, the density of the first block:
                   ptPDens = ptPDens + mL.pDens;
                else
                   ptPDens = ptPDens + mL.sum;
                end

                disp(['For transf. ' parName '=' num2str(thisp) ' Other=' num2str(blockN) '  ,  block sum LL: ' num2str(mL.sum)]);
                ptSLL   = ptSLL + mL.sum; 
                
                if blockN == length(expD{ptN}.other)
                   % NOTE THE PARAM VALUE IS DELIBERATELY NOT STORED HERE.
                   expGrid{ptN,par2grid}.parPDens(2,iGr+grN+1) = ptPDens;
                   expGrid{ptN,par2grid}.parPDens(3,iGr+grN+1) = ptSLL;
                   disp(['param. val.: ' num2str(thisp) '  sum post. dens. for all Others: ' num2str(ptPDens)]);
                end
                
                %%  Now calculate prediction errors that will be used to update
                % pH0, pS0 and d for the next block. NB these are notated D in the
                % 'reporting spokes' but derive from updated beliefs d of the 'hub'
                % for HI:
                post = mL.MDPHI(end).D{1}'; post = post/sum(post); 
                DH0 = sum(post .* midGrid) - pH0;
                % for SI:
                post = mL.MDPSI(end).D{1}'; post = post/sum(post); 
                DS0 = sum(post .* midGrid) - pS0; 
                % d doesn't need this, in this simplest of formulations
                
            end
          
        end
        
        %% ~~~~~~~~~~~~~ Find max and Save as we go along ~~~~~~~~~~~~~~~~ 
        [MAPv, MAPi] = max(expGrid{ptN,par2grid}.parPDens(2,:)); 
        optp = expGrid{ptN,par2grid}.parPDens(1,MAPi);
        if par2grid <= totParN-2
           optPar(ptCnt, par2grid)  = optp;
           refMapMod.mdpStruc.indexP.(fields4MDPs{par2grid}) = optp; 
        elseif par2grid == totParN-1  % this should be othLRTr
           optPar(ptCnt, totParN-1)  = optp;
           refMapMod.mdpStruc.othLRTr = optp;
        else  % this should be pocBTr
           optPar(ptCnt, totParN)  = optp;
           refMapMod.mdpStruc.pocBTr = optp;       
        end
        optPar(ptCnt,totParN+1) = MAPv;   % the fit measures. First is MAP density, then ...
        sumLL = expGrid{ptN,par2grid}.parPDens(3,MAPi); % sum LL, then,
        optPar(ptCnt,totParN+2) = sumLL;
        optPar(ptCnt,totParN+3) = -2*sumLL + 2*par2fitN; % Akaike IC (AIC)
        % BIC := ?2 ln Lopt + dof(params) ln(N/2*pi) = ?2 ln Lopt + dof(params)*(ln N - ln(2*pi))
        optPar(ptCnt,totParN+4) = -2*sumLL + par2fitN*(log(datNPerPt) - log(2*pi)); % Bayesian IC
                                                                       % BIC for small data sample       
        nDone = (ptCnt - 1)*length(parOrdGrid) + parCnt;
        nTot  = length(parOrdGrid)*length(todo);
        disp(['Done ' num2str(nDone) ' out of ' num2str(nTot) ' since ' tStart ]);                                                             
        disp(['Best ' parName ' so far: ' num2str(optp) ' MAP density:' num2str(MAPv)]);
        disp([' All best: ' num2str(optPar(ptCnt,:))]);

        save([fitFileName '.mat'],'expGrid','optPar');
    end
    
    mat2csv2Dfl(optPar,[fitFileName '.csv'],0,1,hd)

end

disp([' attrssri ' wave ' Grid* ' num2str(todo) 'done.'])

return;  % whole function attrssri_Grid10j_vii.




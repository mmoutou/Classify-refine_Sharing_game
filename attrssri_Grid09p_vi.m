function optPar = attrssri_Grid09p_vi( block, wave, startFitFile, grN, grSc, blocksize, codeTest )
%  optPar = attrssribGrid09l_vi( block, wave) :  MAP fitting, with weakly informative priors given by 
% attrssriaGrid09l_vi_**  MAP fitting, with weakly informative priors given by 
%                        scaled beta distributions. ssria means the first wave, and 09l means
%                        that it's the 09 version of the generative function, l for a l set
%                        of parameters to be fitted/considered EXCLUDING POC BIAS and mem) 
%                        Fitting being done with simplest learning from one block to the other,
%                        adding admixture of the posterior to the prior.
%                        Estimate over a grid. vi is the version of the likelihood fn used.
% Modify this at will and just store the version used for a bunch of results in the respective dir.
% To add / remove params to grid over, change:
% 1. par2fitN  2. parOrdGrid 3. fitFileName 3. If need be, if par2grid <= par2fitN-1... statementS
% NB hd and totParN NOT to be changed for this purpose.


% BASIC PIPELINE:
%     . serialDictator* -> formulate and practice MPD based model
%                  and assemble MAP structure / DCM structure for model-fitting to be 
%                  used within EstimParAcIn below (in this fn).
%     . spm_mdp_L_*     -> Likelihood fn 
%     . [experiment]Dat4AcIn*   -> bring exprerimental data to format of generative model
%            e.g. ...\Dropbox\BASOR\AcInSOR\AcInRepeatedDictator\attrssrib\attrSSRIbDat4AcIn09.m
%-->  [ . [expt]Grid[number][let]_[roman]_[pts block}  -> Grid fit INCL learning from one block to the next
%  or { . EstimParAcIn*   -> Use fmincon (or spm_nlsi_Newton...) to fit model structure
%     { . [experiment]Fit[number][letter]_[latin]_[number], e.g. attrssriaFit07a_iv_01 to fit 1st batch of data,
%               with param combination 'a', using likelihood function iv, gen model 07.
%     . mergeExpFits[number][letter]_[number] --> produce nice csv out of multiple fits that will have been 
%               done in parallel, including descriptive data.
fitStr = '09p';

%% Cater for Fixed value parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% column headings for key outputs:
%          1     2      3          4           5        6          7    8     9    10     11    12     13   14    15    16
hd =    {'pH0','pS0','dInitEv','aInitEv','initEvRat','alphaPrec','mem','wH' ,'wS','w0','othLR','pocB','LP','LL','AIC','BIC'};
hd4outp = hd;   % we'll need a backup ;)

parNOTTOGRID = [12, 10]; 
pocBTrRef = log(1.0); 
w0TrRef = 0.0; 

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
try
 load(startFitFile ); % should provide ldat which for 09j should have 29 col and
      % <= 73 rows for wave a, <= 66 for b. It will over-write hd, so fix this below.
 hdLdat = hd;         hd = hd4outp; 
catch
 ldat = [];
end
if strcmp(wave,'a'); totPtN = 73; else totPtN = 74; end
if ~isempty(ldat)
 if size(ldat,1) > totPtN;  error(['too many rows, > totPtN, in ' startFitFile]);      end
 if size(ldat,2) ~= 29; error([startFitFile ' should provide ldat w 29 cols']); end
end

try codeTest; catch codeTest=0; end  % affects grid size, block size etc.
try grN;  catch, grN=12;    end
try grSc; catch, grSc = 5;  end      % will try grids of: par0 + (-grN:grN)*(abs(par0)+1)/grSc ; 
if codeTest
   warning('running with code testing settings; For debug, try grN=1 and grSc either 2 or 200')
end

% NB each pt has 4 partners !! So blocksize = 3 means 12 fits,
try 
    blocksize; 
catch 
    if strcmp(wave,'a')
       blocksize = 18;  
    else 
       blocksize = 16;
    end
end
if codeTest; blocksize = 2; end
todo = (block-1)*blocksize+1 : block*blocksize;
% adjust last block which mops up a bit more:
if strcmp(wave, 'a')
   if (73-todo(end)) < blocksize; todo = [todo,(todo(end)+1):73]; end
elseif strcmp(wave, 'b')
   if (66-todo(end)) < blocksize; todo = [todo,(todo(end)+1):66]; end
else
   error(['wave ' wave '??? ']);
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cwd = pwd;
wrkDri = ['C:\Users\mmpsy\Dropbox\BASOR\AcInSOR\AcInRepeatedDictator\attrssri' wave '\'];
try
   cd(wrkDri);
catch
    try
      wrkDri = ['/home/mmoutou/Dropbox/BASOR/AcInSOR/AcInRepeatedDictator/attrssri' wave '/'];
      cd(wrkDri);
    catch
       wrkDri = ['/home/michael/Dropbox/BASOR/AcInSOR/AcInRepeatedDictator/attrssri' wave '/'];
       cd(wrkDri);       
    end
end

% Loads 'modStruc', which has the key elements of the model structure,
% and expD09 which has the real data, and Pall_template which is a Pall that 
% will be used as a template or scaffold for all the fixed params.
AIDatFile  = ['attrssri' wave 'AcIn09.mat'];            
load(AIDatFile);        % has expD09. Data must be (re-formated) to fit w model Inp, Resp etc.
expDatFile = ['attrssri' wave 'Dat.mat'];               
load(expDatFile);       % This doesn't depend on the generative model.
modStrucFile = 'modStruc09b.mat';              load(modStrucFile);     % This should contain modStruc
datNPerPt = 3 * 4 * 12;  % {data = fairnes prediction, HI, SI ratings} x 4 partners x 12 trials per partner
% also load how the order of wave 2 data corresponds to wave 1:
load('ptIDsAB.mat');

fitFileName = ['exp_' wave '_Grid' fitStr '_' num2str(todo(1)) 'to' num2str(todo(end))];
if grN ~= 12 || grSc ~= 5    % if not the defaults
    fitFileName = [fitFileName '_grSize' num2str(2*grN+1) 'grSc' num2str(grSc)];
end

expGrid = {};
indexPLen = length(fields(modStruc.indexP)); 
par2fitN=indexPLen-1;   % We would like to have the choice to fit any of the params in modStuc.indexP
       % but not nessarily all of them. e.g. length(fields(modStruc.indexP)) = 10
       % but only fit 1st 6, so (indexPLen)+1 has all in indexP plus
       % the 7th param is the leaning rate from partner to partner, othLR
totParN=indexPLen+2; % total number of params in indexP, plus lrnR explored here, AND pocB which is constant.
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
    % in 09p, missing 12 (pocB) and 7 (mem) and 10 (w0=0):
    parOrdGrid = [1 2 8 3 5 11 4 6 9 1 2 8 3 5 6 4 9 1 2 11];  % in order above in hd. Cycle twice!
    if codeTest; parOrdGrid = [8 11 2 11]; end    
    % was: parOrdGrid = [1 2 8 3 5 11 12 4 6 7 9 10 1 2 8 3 5 12 7 6 4 9 1 2 11 12];  % in order above in hd. Cycle twice!
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

           % Params 1 and 2 - used inside MDPs :
           blockN = 1; 
           try
              natPinit.pH0 =   (expD{ptN}.partner(blockN).descr(6)-1)/4;   % rescaled from 1-5 to 0-1 
              natPinit.pS0 =   (expD{ptN}.partner(blockN).descr(7)-1)/4;   % rescaled from 1-5 to 0-1 
           catch  
              warning('Could not set initial values from data - BEWARE!');
              natPinit.pH0 = 0.2;    natPinit.pS0 = 0.2;  
           end

           % Now fill in the rest of the initial conditions for the fit, or 
           % if initial cond. file has been given well, use the values from that.
           if isempty(ldat)    % if no good init. cond. loaded
                % params 3 to 10 - used inside MDPs :
                natPinit.dInitEv = 1;      % 5: try starting from a strongish value!
                natPinit.aInitEv = 1;      % 5: try starting from a strongish value!
                natPinit.initEvRat  =1;    % boring value ...
                natPinit.alphaPrec =0.5; % 0.5: weak value ...
                natPinit.mem  =0.9999;   % no forgetting default
                natPinit.wH  = 6;        % wS and wH added for version attrssriaFit07c_iv_**         
                natPinit.wS  = 6;        % 
                natPinit.w0  = w0TrRef;  % fixed for 09p
           
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
           else    % if good init cond. loaded
% was:            if ~strcmp(hdLdat{18},'othLR'); error('hd{18} ~= othLR'); end
%                 othLRTr = ldat(ptN,18); 
%                 mapMod.mdpStruc.othLRTr = othLRTr ; 
%                 pocBTr = ldat(ptN,19);
%                 mapMod.mdpStruc.pocBTr = pocBTr; 
%                 if ~strcmp(hdLdat{8},'pH0'); error('hd{8} ~= pH0'); end
%                 mapMod.mdpStruc.indexP = modStruc.indexP;  % THIS IS A FILLER CONSTRUCTOR
%                 for k=1:length(spm_vec(mapMod.mdpStruc.indexP))
%                     mapMod.mdpStruc.indexP.(hdLdat{7+k}) = ldat(ptN,7+k);
%                 end
                if ~strcmp(hdLdat{18},'othLR'); error('hd{18} ~= othLR'); end
                if strcmp(wave,'a') % In the first wave, the rows of ldat correspond to ptN,
                            % but in wave 2 they don't. 
                   row4init = ptN;
                else  % we are in wave 2
                   row4init = baIDmap(ptN,2);
                end               
                othLRTr = ldat(row4init,18); 
                mapMod.mdpStruc.othLRTr = othLRTr ; 
                pocBTr = ldat(row4init,19);
                mapMod.mdpStruc.pocBTr = pocBTr;

                if ~strcmp(hdLdat{8},'pH0'); error('hd{8} ~= pH0'); end
                mapMod.mdpStruc.indexP = modStruc.indexP;  % THIS IS A FILLER CONSTRUCTOR
                for k=1:length(spm_vec(mapMod.mdpStruc.indexP))
                    mapMod.mdpStruc.indexP.(hdLdat{7+k}) = ldat(row4init,7+k);
                end

                mapMod.mdpStruc.indexP.w0 = w0TrRef;  % fixed value
           end

           fields4MDPs = fieldnames(mapMod.mdpStruc.indexP); 
           refMapMod = mapMod;     % master copy
        end
        
        %% Iterate over the grid
        for iGr = -grN:grN      
            mapMod = refMapMod; % REM mapMod is max a posteriori focused model whose elements
                                % will be modified acc. to learning, coalitional bias etc. and
                                % then be directly passed to likelihood (incl. MAP density) function.

            %% Parameter value acc. to grid, but before adjustment for which block we're in
            if par2grid <= totParN-2  % i.e. if it's not the learning-between-others learning rate
                refp = refMapMod.mdpStruc.indexP.(fields4MDPs{par2grid}); 
                thisp =  refp + iGr * (abs(refp) + 1)/grSc ; 
                mapMod.mdpStruc.indexP.(fields4MDPs{par2grid}) = thisp; 
                parName = [fields4MDPs{par2grid} 'Tr'];
            elseif par2grid == totParN-1     % i.e. if stepping the others-lr
                refp = othLRTr;
                thisp =  refp + iGr * (abs(refp) + 1)/grSc;                 
                mapMod.mdpStruc.othLRTr = thisp; 
                parName = 'othLRTr';
            else   % Shouldn't get to this!
                disp('  OUCH '); 
                error(['Why did we get to par2grid=' num2str(par2grid)]);
                refp = pocBTr;
                thisp =  refp + iGr * (abs(refp) + 1)/grSc ; 
                mapMod.mdpStruc.pocBTr = thisp; 
                parName = 'pocBTr'; 
            end
            
            %% ENSURE PARAM VALUES NOT TO BE EXPLORED ARE FIXED
            % pocB set to neutral
            pocBTr = log(1.00);                             % fixed unbiased value
            mapMod.mdpStruc.pocBTr = pocBTr;                % ... to real line ... 
            mapMod.mdpStruc.pocBTr = pocBTr;                % ... to real line ... 

            % mem set to no-forgetting, logit(0.9999)
            mapMod.mdpStruc.indexP.mem = 9.210240366975960; 
            % w0 set to 0 (again - redundant?)
            mapMod.mdpStruc.indexP.w0 = w0TrRef;
            %%
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
            aInitEv = exp(mapMod.mdpStruc.indexP.aInitEv); 
            othLR = 1/(1+exp(-mapMod.mdpStruc.othLRTr));

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
                    aInitEv = aInitEv + lambda;  
                    if dInitEv < 0.001; dInitEv = 0.001; warning('dInitEv set to 0.001 :( '); end
                    if aInitEv < 0.001; aInitEv = 0.001; warning('aInitEv set to 0.001 :( '); end
                    mapMod.mdpStruc.indexP.dInitEv = log(dInitEv); 
                    mapMod.mdpStruc.indexP.aInitEv = log(aInitEv);                     
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
                mL = spm_mdp_L_vi(mapMod.mdpStruc.indexP,mapMod.mdpStruc,...
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
                
                %%  Now calculate PEs (here, belief shifts) that will be used to update
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
        else  % this should be NOT pocBTr
           error('Why did we get here??');
           optPar(ptCnt, totParN)  = optp;
           refMapMod.mdpStruc.pocBTr = optp;       
        end
        refMapMod.mdpStruc.pocBTr = mapMod.mdpStruc.pocBTr;  % pocB is constant in 09k.
        optPar(ptCnt, 12) = refMapMod.mdpStruc.pocBTr;       % (some redundant code here for recording same ...)
        
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
        
        save([fitFileName '.mat'],'expGrid','optPar','hd');
    end

    mat2csv2Dfl(optPar,[fitFileName '.csv'],0,1,hd)

end

disp(' '); 
disp(['Grid for  attrssri ' fitFileName ' done.'])
t = clock; tEnd = num2str(t(2:5)); 
tStart
tEnd
save([fitFileName '_done.mat'],'tStart','tEnd','parOrdGrid','grSc','grN');

return;  % whole function.





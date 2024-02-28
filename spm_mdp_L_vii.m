function mL = spm_mdp_L_vii(P,modStruc,Inp,Resp,details)
% minus-log-likelihood function, developped from spm_mdp_L to go with serialDictator10a - roman num index _i, _lxix etc.
% FORMAT mL = spm_mdp_L_vii(P,modStruc,Inp,Resp)
% P         - parameter VECTOR or structure, IN -INF TO INF (transformed not native) space
%             This fun. must be able to use with fminunc etc. which expect P to be vector.
% modStruc  - generative model (or sufficient ingredients thereof)
%             It is strongly advised for modStruc to include modStruc.indexP as a STRUCTURE according to 
%             which we will cast a vectorial P.
%           - Also can have modStruc.priPar, with the (macro)params of the prior distros 
%             in rows of alpha, beta, lo, hi for use by dbetasc .
% Inp aka U - inputs to gen. mod.   ('obsevations': stimuli, outcomes, where we are ...)
% Resp aka Y - observed outputs       (responses / actions performed)
% Resp, or both Resp and Inp, can be empty matrices if we want to run in generative mode. 
%
% BASIC PIPELINE:
%     . serialDictator* -> formulate and practice MPD based model
%                  and assemble MAP structure / DCM structure for model-fitting to be 
%                  used within EstimParAcIn below (in this fn).
% ->  . spm_mdp_L_*     -> Likelihood fn 
%     . [experiment]Dat4AcIn*   -> bring exprerimental data to format of generative model
%     . EstimParAcIn*   -> Use fmincon (or spm_nlsi_Newton...) to fit model structure
%     . [experiment]Fit[number][letter]_[latin]_[number], e.g. attrssriaFit07a_iv_01 to fit 1st batch of data,
%               with param combination 'a', using likelihood function iv, gen model 07.
%     . mergeExpFits[number][letter]_[number] --> produce nice csv out of multiple fits that will have been 
%               done in parallel, including descriptive data.
%
% This dovetails with serialDictator09b
% This function finds the key ingredient(s) of the generative model  in modStruc
% and a given set of (transformed to -inf to inf) parameter values in P and 
% values, after adding in the observations and actions on each trial
% from (real or simulated) participant data. It then sums the
% (log-)probabilities (log-likelihood) of the participant's actions under the model when it
% includes that set of parameter values. 
%
% Fitting routines can use this function to fit model params to data. 
%
% Demo/test with: 
%                       [ HI SI]         [pH0,pS0,initEv,uPrec,initEvRat,wH,wS, mem, desBias,desCorr] 
% Demo: close all; pOth=[0.9,0.4]; pSel= [0.1,0.3, 2,   0.1,    1.5  ,   6, 6,  0.99,  0.0,     2]; 
% close all; pOth=[0.9,0.4]; pSel= [0.1,0.33, 1,  1/5,  1,        6, 6, 0.99,  0.0,   5]; [MDPs,modStruc,Inp,Resp] = serialDictator07(pSel,pOth,1); 
% modStruc.priPar=[[1.01 ,1.01,1.01, 1.01, 1.01 ,10]; [1.01 ,1.01, 2,  2, 2, 10]; [ 0, 0, 0, 0, 0, -46]; [ 1, 1, 100, 100, 100, 46]];   
% P = spm_vec(modStruc.indexP) ;    test = spm_mdp_L_vi(P,modStruc,Inp,Resp,1); disp(test);  P(5) = P(5)/10; test = spm_mdp_L_iv(P,modStruc,Inp,Resp,1); disp(test); 
%__________________________________________________________________________
debugging = 0;   % will report total sum LL + logPrior if 0, or just the Lattribution if 2,
                 % or just the Lprediction if 1
try 
    details;     % details > 0 enriches mL, >= 2 plots stuff.
catch
    details=0;                     % default if actions given is to just output the log lik
    if isempty(Resp); details=1; end  % ... but actions not given, output lots of synthesized stuff.
end   


Pall   = modStruc.allP;    % comprehensive list of params to modify and 
                           % use to change d, A, C and other maps
tiny = Pall.tiny;
trialN = Pall.trialN;      resNRepo = Pall.resNRepo; 
hmmHub = modStruc.hmmHub;    
% The following lines were needed if there was going to be a learning:
%     Make sure that there is no A field if we are running in generative mode
%     and inputs (pos/unfair returns etc) have been provided: 
%     if isfield(hmmHub,'A') && ~isempty(Inp); hmmHub = rmfield(hmmHub,'A'); end
    
mdpFair = modStruc.mdpFair;
mdpHI = modStruc.mdpHI;      mdpSI = modStruc.mdpSI;    

%% Bring parameters to native space, i.e. out of log- , logit- etc. space  
% to make suitable to pass to the model to compute the log-likelihood
%--------------------------------------------------------------------------

% REM P AND indexP have ONLY the 'index' params, e.g. to be fitted! 
%              Not the ones needed for additional specification / fixed ones! 
%              To see all the params, incld he fixed ones, list modStruc.allP .
% First, bring P from its inputted form, such as a vector, to the form that 
% the param structure modStruc.indexP also has :
if ~isstruct(P); P = spm_unvec(P,modStruc.indexP); end

% default is not to adjust the various maps unless dictated by P.
% Here record only the maps that may be meaningfully adjusted.
% Learning 'hub' :
dHubAdj=0;      CHubAdj = 0;    %  aHubAdj=0;
% Rep.orting 'spokes':
ARepAdj=0;      CRepAdj=0;      pFHSAdj = 0;

% Complicated if statement to transform inputted parameters to native 
% space and replace the relevant values in Pall
field = fieldnames(P);   
Plen = length(field);   

for i = 1:Plen
    % first, log-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if strcmp(field{i},'dInitEv')
         Pall.dInitEv = exp(P.dInitEv);       dHubAdj=1;
%      elseif strcmp(field{i},'aInitEv')
%          Pall.aInitEv = exp(P.aInitEv);       aHubAdj=1;
     elseif strcmp(field{i},'initEvRat')
        Pall.initEvRat = exp(P.initEvRat);    dHubAdj=1;  % aHubAdj=1;
     elseif strcmp(field{i},'alphaPrec')
        Pall.alphaPrec = exp(P.alphaPrec);   
        hmmHub.alpha = Pall.alphaPrec;      
        mdpHI.alpha = Pall.alphaPrec;   
        mdpSI.alpha = Pall.alphaPrec; 
        mdpFair.alpha = Pall.alphaPrec; 
     elseif strcmp(field{i},'wH')
        % Here we check before changing, bec. if we don't 
        % need to change we can save an expensive function call later.
        exp_wH = exp(P.wH);
        if abs(Pall.wH - exp_wH) >  tiny
           Pall.wH = exp_wH;    
        end
    elseif strcmp(field{i},'wS')
        exp_wS = exp(P.wS);
        if abs(Pall.wS - exp_wS) >  tiny
           Pall.wS = exp_wS;   
        end    
    elseif strcmp(field{i},'Ucor')
        Pall.Ucor = exp(P.Ucor);            ARepAdj=1; 
    % logit-transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'pH0')
        Pall.pH0 = 1/(1+exp(-P.pH0));       dHubAdj=1;
    elseif strcmp(field{i},'pS0')
        Pall.pS0 = 1/(1+exp(-P.pS0));       dHubAdj=1;
    elseif strcmp(field{i},'lrnR')  %  v. simple - change directly:
        Pall.lrnR = 1/(1+exp(-P.lrnR));  
        hmmHub.eta = Pall.lrnR; 
    elseif strcmp(field{i},'mem')  %  v. simple - change directly:
        Pall.mem = 1/(1+exp(-P.mem));  
        hmmHub.omega = Pall.mem; 
    % and scaled-logit transformed: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'desBias')
        Pall.desBias= -1+ 2/(1+exp(-P.desBias));   ARepAdj=1;  
    % Untransformed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(field{i},'w0')
        if abs(Pall.w0 - P.w0) >  tiny
         Pall.w0 = P.w0;        
        end
    elseif strcmp(field{i},'desCorr')
         Pall.desCorr = P.desCorr;   CRepAdj=1;  
    elseif strcmp(field{i},'desFair')
         Pall.desFair = P.desFair;   CHubAdj=1;   
    else
       error([field{i} ' not catered for :(']);
   end
end

%% Prior density over now adjusted parameters:
%  NB if we want to ignore log pri over specific param theta, simply
%  set them sampled from betasc(theta,1,1,-huge,huge)
try priPar = modStruc.priPar; catch, priPar = []; end
if ~isempty(priPar)  % will also work if modStruc.priPar exists and was empty.
    natPvec = zeros(1,Plen);   % for P in native space
    for i = 1:Plen
        natPvec(i) = Pall.(field{i});
    end
    logPri = log( dbetasc(natPvec,priPar(1,:),priPar(2,:),priPar(3,:), priPar(4,:)));  
else
    logPri = 0;
end


%% ** Assemble and run the generative model with correct stimuli and returns  **

%% Extract the single-trial MDPs from modStruc, components of the overall model
%  that we will need to put together. Also extract other auxiliaries.

% Adjust component MDP maps acc to parameters as necessary  -------------
%% for Hub:
if   CHubAdj; hmmHub.C{1}  = CHubSerDict(Pall);     end 

if  dHubAdj 
    d0 = d0HubSerDict(Pall);   % Will need copy before it's updated!
    hmmHub.d     =  d0; 
end

%% For spokes (report MDPs)
if  ARepAdj
    mdpHI.A{2} = ARepClassifSerDict(Pall);   % The 'Classif' here has larger resoln. than in _vi
    mdpSI.A{2} = mdpHI.A{2};  
    mdpFair.A{2} = mdpHI.A{2}; 
end
if CRepAdj
    mdpHI.C{2} =  CRepClassifRepDict(Pall); 
    mdpSI.C{2} =  mdpHI.C{2};
    mdpFair.C{2} =  mdpHI.C{2};
end

%%  1. put together and run 'Hub', or Learning MDP ------------------------------------

% Assemble trials and add observations and actions ------------------------------------
[hmmHub(1:trialN)] = deal(hmmHub);          % Create MDP with specified number of trials

% Cater for deliberately missing function arguments, when fn. used in 
% generative mode : 
if ~isempty(Inp)
   if isempty(Resp)    % if the actions have not been provided, they will be estimated ...
       % ... but unless we set the entries below to 0, the corresponding observations,
       % which are contingent on actions, will be 'stuck' in those in Inp. However they 
       % **should** be stuck to those in Inp if Resp has been provided :)
       for k=1:trialN;  Inp.hub{k} = 0;  end
   end
   [hmmHub.o]      = deal(Inp.hub{1:trialN});     % Add observations in each trial
end

%--------------------------------------------------------------------------
%  Run model for hub, so as to map beliefs about each available state
%  in each trial, given parameters and stimuli :  
hmmHub   = spm_MDP_VB_XI(hmmHub); % run model with given parameter values
%--------------------------------------------------------------------------

% Store true generated states and flag appropriately if in generative mode:
if isempty(Resp)
    Resp.sHub = {hmmHub.s};     
    genMode = 1;
else
    genMode = 0;
end   

%%  Reporting MDPs ---------------------------------------------------------------------
% Observations / stimuli in this model include observations of one's
% own responses, as well as fair or unfair returns. Only the third row
% of observations of the 'learning MDP' contains the fair/unfair ret.,
% and this does not depend on pts. beliefs and actions, hopefully!
% Should always go neutral - fair_or_unfair 

% Accumulate log-likelihood for Reporting MDPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REM: MDP.P(M1,...,MF,T)    - probability of emitting action M1,.. over time

LAttr  = 0;  LFair=0; % initialise summed (log) probability of reporting MDP actions

for tr = 1:trialN
    
  % Copy from updated-parameter element prototypes:
  mdpH = mdpHI;          mdpS = mdpSI;      mdpF = mdpFair; 
  
  % For the d map of the Fairness reporting MDP, we must use the posterior a 
  % of the PREVIOUS trial, or the very starting one:
  if tr > 1
     % mdpF.d = hub2spoke_dFair(hmmHub(tr-1).a{1}, hmmHub(tr-1).d, Pall);
     mdpF.D = hub2spoke_DFair_vii(hmmHub(tr-1).d, Pall);
  else
     mdpF.D = hub2spoke_DFair_vii(d0, Pall);  % use copies from before hmmHub model ran at all!
  end
  
  % Set true states for the attribution MDPs. REM columns are within-trial timepoints, 
  % rows are state factors. So top row is HItrue, etc.
  mdpH.s(1,1)   = hmmHub(tr).s(1,1);       mdpS.s(1,1)   = hmmHub(tr).s(2,1);
  mdpH.s(2,1)   = resNRepo+1;              mdpS.s(2,1)   = resNRepo+1;
  
  % Set the attribution d maps according to the updated hubHMM. Here, we need 
  %    to know the posterior d for the hub MDP, to see what attributes the new
  %    beliefs about d amount to. To derive what the 'hi' and 'lo' states
  %    now amount to, we'll  use the portion of aHub that maps 
  %    from hi-lo HI and SI states to returns.
  [mdpH, mdpS] = hub2spoke_DAttr_vii(mdpH, mdpS, hmmHub(tr), Pall); 
  % Test with e.g. squeeze(mdpHI.A{2}(:,3,:)) to see correctness levels for all Hreport when Hbelieved is 3

  % Set the crucial actions: 
  if ~genMode
      mdpF.u = Resp.fairRep{tr};
      mdpH.u = Resp.Hattr{tr}; 
      mdpS.u = Resp.Sattr{tr}; 
  end
  
  %% ~~~~~~~~~~~~~~~ Solving and Storing the Attibuting MDPs ~~~~~~~~~~~~~~~~~~~~~
  
  % Now solve and store:
  mdpH  = spm_MDP_VB_XI(mdpH);             mdpS  = spm_MDP_VB_XI(mdpS); 
  mdpF  = spm_MDP_VB_XI(mdpF); 
  
  if genMode || details % i.e., if we are in generative mode or user wants gory details
    if genMode % test again, lol
      Resp.fairRep{tr}  = mdpF.u;    
      Resp.Hattr{tr} = mdpH.u;    
      Resp.Sattr{tr} = mdpS.u;        %#ok<*SAGROW>
    end
    if tr == 1
      MDPH(1:trialN) = deal(mdpH);  % }  These store first trial and do memory prealloc.
      MDPS(1:trialN) = deal(mdpS);  % }
      MDPF(1:trialN) = deal(mdpF);  % }
    else
      MDPH(tr) = mdpH;     MDPS(tr) = mdpS;    MDPF(tr) = mdpF; %#ok<*SAGROW>
    end
  end
  
  % clear mdpH mdpS mdpF; % inefficient but safe ...

  % Get probability of  action to be taken for each reporting MDP:
  l1 = log( mdpF.P(1, Resp.fairRep{tr}(2,1)) + 1e-16);   
  l2 = log( mdpH.P(1, Resp.Hattr{tr}(2,1)) + 1e-16); 
  l3 = log( mdpS.P(1, Resp.Sattr{tr}(2,1)) + 1e-16); 

  LAttr = LAttr + l2 + l3;
  LFair = LFair + l1; 
    
  if details  % Store progress in detail for debugging etc. ~~~~~~~~~~
        mL.lFair(tr) = l1;     mL.lHI(tr) = l2;      mL.lSI(tr) = l3;                  
  end % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
end

if details  % Detailed output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % record stuff
    mL.sum = LFair + LAttr;                 % likelihood   
    mL.pDens = mL.sum + sum(logPri);      %  posterior density, for MAP estimation 
    mL.Lfair=LFair;            mL.LAttr=LAttr;                
    mL.hubHMM=hmmHub;   mL.MDPFair = MDPF;       
    mL.MDPHI = MDPH;    mL.MDPSI = MDPS; 
    mL.P = P;                  mL.Pall = Pall; 
    mL.Inp.hub = {hmmHub.o};   mL.Resp = Resp;   % Copies handy if we ran in generative mode.
    mL.logPri = logPri;        % log priors over parameters.
    mL.d0Hub = d0;
    
    % plot stuff if asked:
    if details >= 2
      
      spm_figure('GetWin','Attrib. Harm Intent, trial-by-trial'); clf  
      spm_MDP_VB_game(MDPH);
      spm_figure('GetWin','Attrib. Selfish Intent, trial-by-trial'); clf  
      spm_MDP_VB_game(MDPS);
	  spm_figure('GetWin','Fairness reporting: All, trial-by-trial'); 
	  spm_MDP_VB_game(MDPF);
    
      % aHubProbDisp displays how the likelihoods in a are learnt over trials,
      % but here, in version vii, there is no a learning
      % aHub = mdp2arr(hmmHub,'a');  aHub = aHub{1}; 
      % aHubProbDisp
         
    disp(' ');
    disp('~~~~~~~~~ Key to muliple trial plots: ~~~~~~~~~~~')
    disp('Topmost  : Circles: true init. state along each factor, from mdp.s. Colours are codes:')
    disp('  Modulo numeric <-> {''r.'',''g.'',''b.'',''c.'',''m.'',''k.''}; ')
    disp('  Rows of circles corresp. to factors, first factor TOP-MOST, ')
    disp('  e.g. if initial fairnessLevelReportState=5 -> magenta')
    disp('2nd down : Observations / outcomes for each factor, from mdp.o ')
    disp('  TOP-MOST is LAST factor here :(, fairness of split, 3=null/not present in timestep, 1=UNFAIR')
     % disp('           ')
    disp('3rd down : ERP ;     4th down : DA ')
    disp('5th down : Post. beliefs about states, (D / updated d). y axes bizarelly seem to go HIGHER DOWN! ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
    end
    
    clear('hubHMM','MDPHI','MDPSI'); 
% Else the minus sum log likelihoods for different types of fit .... ~~~~~~~
elseif debugging == 0
    mL = -(LFair + LAttr + sum(logPri));
elseif debugging == 1                   % to fit just fairness report learning
    mL = -(LFair + sum(logPri));
elseif debugging == 2
    mL = -(LAttr + sum(logPri));
else
    error(['Debugging option ' num2str(debugging) ' not catered for']);
end


% clear('hubHMM','mdpHI','mdpSI');    

return;

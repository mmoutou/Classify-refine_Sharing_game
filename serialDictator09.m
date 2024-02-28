function [MDPs, modelStructure, Inp, Resp] = serialDictator09(selfp, otherp, toPlot) 
% serialDictator09 - reporting chance of fair return, updating
%        beliefs based on actual return, AND attribution reporting. 
%        Version to use good-bad classification plus A learning, after
%        Pike 'exponentially reducing learning rate'.  
%        09 is first with 3 independent reporting MDPs, 2 x attrib and 1 x fairness reporting
%        MDPs draw inputs from the learning-MDP, (up to 06 had 1 attribution MDP).
%        This factorises the solution and means that there is no feedback from
%        'HI/SI attributions I will/have declared' to 'underlying beliefs about HI/SI'.
%        toPlot =0: no plots; 1: behaviour & key beliefs; 2: key neural simulations.
%        version 06 explored 6 levels of reporting, 'boring' wH ans wS,
%        'boring' desBias (could be used in AClassifAttrSerDict), but added dInitRat
%                          [ HI SI]         [pH0,pS0,initEv,uPrec,dInRat,wH,wS, w0,lrnR,desCorr, desBias] 
% Demo: close all; otherp=[0.5,0.6]; selfp= [0.4,0.7, 1,    1/5,    1,    6,  4,0.5,0.1,  0.5,   0]; [MDPs,modStruc,Inp,Resp] = serialDictator09(selfp,otherp,1);  
% Uses hidden states causing (factorised) outcomes. Here the factorisation is
% explicit, enabling us to model multiple modalities (outcome factors) and
% distinct hidden causes of observation (hidden state factors like what and
% where). 
%______________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging, Karl Friston et al

testing = 1; % 0 -> use time-based random gen ...
seed = 661;    % Only used if testing~=0

% set up and preliminaries
%==========================================================================
if ~testing
    rng('shuffle'); 
else
    rng(seed, 'twister');
    warning('serialDictator08 running with testing/debug settings incl rng(seed,''twister''');
end % for reproducible sequences of random numbers; 'shuffle' bases on time.

%% Menu - like items -----------------------------------------------------
% Number of trials (maybe in future per each Other encountered)
trialN = 12; % 36;
try selfp;  catch, selfp  = []; end
try otherp; catch, otherp = []; end
if isempty(otherp);  otherp = [0.9, 0.4];   end  % [truHI truSI ]
trueHI = 1+1*(otherp(1) > 0.5);  % call these, a little arbitrarily, 'high' and 'low'
trueSI = 1+1*(otherp(2) > 0.5);  % or 1 and 2, suitable for indexing, for s states, etc.

% Initial binary alternatives -  
wBin = logit(0.25);
HIInit = [invlogit(wBin) (1-invlogit(wBin)) ];    % differentiate if need ease ...
SIInit = HIInit;                                  % ... of checking code.

% True 'alternatives'; or pseudo-alternatives, as true state is fixed to one of them:
if trueHI == 1; HIVal = [otherp(1) 0.75];  else HIVal = [0.25 otherp(1)];  end
if trueSI == 1; SIVal = [otherp(2) 0.75];  else SIVal = [0.25 otherp(2)];  end

%                               1   2   3     4      5       6  7   8    9    10      11
%                             [pH0,pS0,initEv,uPrec,dInitRat,wH,wS,w0, lrnR, desCorr,desBias] 
if isempty(selfp);   selfp  = [0.1,0.3, 2,   0.1,    1,      6, 4, 0.5, 0.1, 1,     0.0];   end  
                                 % REM Joe's reversal paper had pH0,pS0,uH,uS, wH,wS,w0,uPi,eta_dg
                                 % uPrec (cf. uPol) is 1/alphaPrec, the expected prior prec.
pH0 = selfp(1); pS0 = selfp(2); initEv = selfp(3);   % give error if this fails :)
try 
    dInitRat = selfp(5);
    wH=selfp(6); wS=selfp(7); w0 = selfp(8); lrnR= selfp(9); 
    %    bias from -1 to +1, to shift pmf; ...] 
    alphaPrec = 1/selfp(4); % selfp(4)=uPrec (cf. uPol) is 1/alphaPrec, the expected prior prec.
    desCorr = selfp(10); 
catch
    wH=6; wS=4; w0=0.5; lrnR=0.1;  alphaPrec = 10; dInitRat = 1;  desCorr = 1;
end
Ucor = 0.2; % Noise param for noisyBino over wrong ... correct. 0.2 -> care, 2 -> don't care
            % See A3 and A4AttrSerDict - How well HI, SI is to be reported.
Upop = 1.0; % see dClassifSerDict -- param for belief of SD of population, reflected  
            %     in bluntness of d. So set to 1 for representative mean, 0.2 for representative mode (?).
desFair=1;  % scaling of desirability to get fair return.
%warning(['desFair set to ' num2str(desFair)]);

% toPlot = 4 to display neuronal phase precesion / theta-gamma etc.
%       or 1-3 for more basic plots.
try toPlot; catch, toPlot =0; end
% neuroFami = 1 to display familiarity / MMN / other_intent learning - NOT READY.
neuroFamil =0; %try, neuroFamil; catch, neuroFamil =0; end

%% Model specifiation: dimensions of matrices and descriptions of said dimensions/factors
Tsteps1 = 2;        % Each stage of learning-MDP, hmmHub will have 2 Tsteps.
                    %  initial state -> receive-return state. 
Tsteps2 = 2;        % Each stage of reporting MDP, will only have just 2 Tsteps.
                    %  initial state -> choose level-to-report 
stateFNHub = 2;       % Number of state factors for hub HMM == learning-MDP: trueHI, trueSI    
stateFNRep = 2;    % N. of state factors for reporting-MDP: beliefState, reportState
% N. of (pseudo) observations of whether we reported well or not. Here using 3 levels 
% of correctness of reporting, e.g. quite correct=3, not too bad, ...
% ...at most 1 off in each dimension =2, wrong ie the rest =1
corLevN = 3;   %  corLevN+1 is the 'indifferent/neutral'
 
resNHub  = 2;   % resoln. level for classifying partner, starting w. e.g.
                % unfair (0.1 return) fair (0.45) altruist (0.8)
resNRepo = 6;   % Likert-like resolution for reporting fairness or attributions - resolN intervals, w. 
                % boundaries (1:(resolN-1))/resolN and midpoints 1/(2*resolN)+(0:(resolN-1))/resolN
             
              
midGridRepo = 1/(2*resNRepo)+(0:(resNRepo-1))/resNRepo;  % grid of midpoints at specified resolution
ubGridRepo  = (1:resNRepo)/resNRepo;                   % grid of upper bounds of same.
resolRepo   = 1/resNRepo; 

% Generative model parameters to implement Other's policy:
desBias = 0;  % Don't use this but w0 to account for overall biases if poss.
              % Can be used (is coded into) in AClassifAttrSerDict
% midPFair is a resoln x resoln grid of 'choice A vs B' probabilities depending 
% on [wH,wS,w0]
midPFair = midpfair([wH,wS,w0],resNHub);
Ioff = 0.5;  % I always use the same offset for the logistic below ...
corePFair = 1/(1+exp(w0+(HIVal(trueHI)-Ioff)*wH + (SIVal(trueSI)-Ioff)*wS)); 
[~, corePBin] = min(abs(midGridRepo - corePFair));  % Bin where the intended generative
                                                    % fair-share proportion resides.

% pFHS provides all the fair return probabilities for all the combinations of HI and SI
% states that the participant considers if HI and SI can be either high or low, and 
% these low and high- attribute states can only take values from attrV below. iLo and iHi
% are the indices over attrV for each row of pFHS, for the HI and SI responsible acc. to attrV
if resNRepo == 6
    [pFHS, attrV, iLo, iHi, iHS] = pFairHiLoHISI6([w0,wH,wS]);  
else
    error('resNRepo~=6 so pFairHiLoHISI6 not appropriate');
end
% aInit   = initEv; % effective initial evidence for the a map (prior on A ...) 
                 % Simplifying assumption that ppl who are sure about what the state
                 % of the world is will be equally sure about what outcomes to 
                 % expect from states of the world. But to make into free par in due!

% Make a grid of index probabilities to aid noisify A2 (attributing) map. Each
% col of pCorLev is the A map for different correctness levels (of a certain report). I.e.
% if I really really want to report HI=2, then col 2 would be
% noisyBino(corGrid(corLevN),0.2,corLevN) . If I thought that reporting HI=2 is OKish, 
% then that col in A would be noisyBino(corGrid(corLevN-1),0.2,corLevN)
% Desirability shifts will then modify these, so if des > 0 then the HI reports which are
%              lesser than the underlying belief about HI will be shifted towards correct, up,
%              and the >= towards wrong. Vice versa if des < 0. 
if corLevN == 3
    corGrid = [0.2, 0.5, 0.8];
else
    corGrid = 1/(2*corLevN)+(0:(corLevN-1))/corLevN;
end
%  end block about wrong - OKish - correct - Neutral stuff 

%% re-record for output the params ----------------------------------------------
% First, the ones that we will want fitted to data - see
%  Likelihood function e.g. spm_mdp_L_*.m , where indexP stands for 
% structure of index params, as transformed to the real line. 
% Also, a more extended structure with all the native-space params that
% that may be needed to construct d0, AHub, A2 etc in spm_mdp_L*.m 
p4fit.pH0   = pH0;              p4fit.pS0 = pS0;  
p4fit.initEv = initEv;          p4fit.dInitRat = dInitRat;
p4fit.alphaPrec = alphaPrec; 
p4fit.w0 = w0;                  % p4fit.aInit = aInit;  % no - try setting all from initEv.
p4fit.lrnR = lrnR;  
modelStructure.indexP = nat2tr_mdp_L_ii(p4fit); % indexP in transformed space.
% So if we assume to-be-fitted to be as above,
allP = tr2nat_mdp_L_ii(modelStructure.indexP);  % i.e. back to native, for funsies, how inefficient is this way of doing it ha ha ;)
allP.wH = wH;                allP.wS = wS; 
allP.desBias = desBias;   
allP.desCorr = desCorr;      allP.desFair = desFair;
allP.resNRepo = resNRepo;    allP.resNHub = resNHub;   
% No. of levels of desirability for unfair, fair, ... indifferent/neutral :
allP.retLevN = 3;
allP.Ucor = Ucor;            allP.Upop = Upop; 
allP.corLevN = corLevN; 
allP.corGrid = corGrid;
allP.trialN = trialN;    allP.Tsteps1 = Tsteps1;  allP.Tsteps2 = Tsteps2;
allP.tiny = 1e-16;        % probabilities of this order are zero for all intends and purposes.
allP.noiseFloor = 0.001;  % General purpose 'small compared to other sources of noise'

% modelStructure.d0 = d0;  % DONE BELOW!
modelStructure.allP = allP;  % Store for output!

%%  ############################################################
%%  ######             Learning :   Hub MDP            #########
%%  ######    observe split & update HI & SI beliefs.  #########
%%  ############################################################

%%  ?? V map [ allowable policies (of depth T, so T-1 rows) ] here!
%   For the Hub, there are no actions to be taken - or one inconsequential one.
%   See function spm_MDP_get_T  in file spm_MDP_VB_XI.m  l. 1323
%   ----------------  so no V or U here ! ----------------------------
%   In general, policies are sequences of actions
%   (with an action for each hidden state factor, enumerated in the 
%    last column. T-1 rows bec. actions not taken in last state.)
%  So of the form V(timestep, action-sequence, state-factor) I think
%  REM here stateFNHub = 2;  % No. of state factors: trueHI, trueSI                       
%--------------------------------------------------------------------------

%% %% The key dynamics: B map, i.e. Transitions conditioned on action %% %
% controlled transitions: B{u}
%--------------------------------------------------------------------------
% We specify the probabilistic transitions of hidden states for each
% state (as opposed to outcome) factor. Hence actions are themselves
% factorised, state factor F being subect to actions MF :
%  REM from spm_MDP_VB_X: MDP.B{F}(NF,NF,MF)  - transitions among states 
%       under MF control states. So each B{k} refers to *state factor* k,
%       which has Nk states; There must be a grand total of Mk 
%       actions / control states for these, the 'pages' of B{k}, not
%       all of which might be available under V for the current time-step.
%--------------------------------------------------------------------------
% B1{1} : transitions over factor 1, trueHI : There is only the
%        the action 'accept fate', the pt can't change this.
%        trueHI is considered to have resNHub poss. levels: 
B1{1} = zeros(resNHub,resNHub,1);
B1{1}(:,:,1)= eye(resNHub);   
% Similarly for trueSI:
B1{2} = B1{1};

%% %%%% Learning MDP outcome probabilities: A (likelihood map) %%%%% %
%--------------------------------------------------------------------%
% REM from spm_MDP_VB_X: MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes 
%     given hidden states. G must stand for the *outcome factors*, whereas
%     F must stand for the *state factors*, N1 being the number of states
%     in state factor 1 etc. 
% Here specify the probabilistic mapping from hidden states
% to observations (a.k.a. outcomes):
% AHub{1} : Observe fairness level reported
% AHub{2} : Actual unfair or fair ret.; depends on Stage, trueHI and trueSI.
% Their form is A{G}(observation, trueHI, trueSI, level-to-report)
% See Emotion_learning_model... A map as example. 
%--------------------------------------------------------------------%

% A : Actual fair ret. depends only on trueHI and trueSI
%                 returns     trueHI    trueSI   
% AHub = zeros(retLevN-1,    resNHub,  resNHub); 
aHub = AHubClassSerDict(allP);                      % initial generative ...
aHub(1:2,:,:) = aHub(1:2,:,:) * allP.initEv;        %          ... model likelihoods
AHub = AHubClassSerDict(allP, HIVal,  SIVal);       % generative process likelihoods


%% Learning MDP priors: (utility) C1 aka goal map ---------------------------------------
%--------------------------------------------------------------------------
% Specify the prior preferences in terms of log probabilities over outcomes. 
%    REM from spm_MDP_VB_X: 
%    MDP.C{G}(O,T) - (log) prior preferences for outcomes (modality G), I
%    believe at each of T timesteps or moves ??
%    C factors indexed by G must corresp. to A.
%--------------------------------------------------------------------------
     
% Outcome is desirability of fair etc. split.
%   allP must include, for use here, just desFair
C1  = CHubSerDict(allP); 

%% Prior beliefs about initial states: d map ------------------------------
% now specify prior beliefs about initial states, in terms of counts. Here
% the hidden states are factorised into 
% MAKE SURE THEY ARE COL VECTORS AS IN E.G. d{2} = [2 2]';
% state factors: stage, trueHI, trueSI, Level2Report 
%--------------------------------------------------------------------------

%  d  has 2 factors, beliefs about 'true HI' state & 'true SI'
%  allP must include, for use here: pH0, pS0, initEv, dInitRat, resNHub. 
%     Optional Upop (default=1), the dispersion param. for the population distro.
%  REM below includes r = allP.dInitRat;  dInitH = allP.initEv * r / (1+r);
%                     dInitS = allP.initEv / (1+r); 
d0 = dHubClassSerDict(allP); 

%%  Hub MDP Structure - to be used to generate arrays for multiple trials
%==========================================================================

%% Generative and process model structure for hub -  -  -  -  -  -
% Note conc. parameters are counts not probs. Specifying the mdp.a here 
% is sufficient for the a map to undergo learning when the model is 'solved' below.

hub_a0{1} = aHub;   % hub_a0 will serve also as a working and storage copy
                    % intialise! Best to keep cell array notation even for 1 element.
hmmHub.a = hub_a0; 
hmmHub.eta = lrnR;    % This is 'learning rate' (learning between trials, I hope).

% The true generative process:
hmmHub.A = AHub; 
hmmHub.B = B1;                  % transition probabilities, allready cell array
hmmHub.C{1} = C1;               % preferred outcomes
hmmHub.d = d0;                  % prior over initial states
hmmHub.s = [trueHI, trueSI]';   % true initial state. [resNHub, resNHub] would mean
                                % highest true HI and SI
hmmHub.alpha = alphaPrec;
hmmHub.T = Tsteps1 - 1;         warning('Why hmmHub.T = Tsteps1-1 works??'); 
hmmHub.tau   = 12;
%  fix the names
hmmHub.Aname = {'returnFairness'};
hmmHub.Bname = {'trueHI','trueSI'};
% The following line doesn't work properly - there may be a bug in spm_MDP_VB_trial.m
hubModalityNames =  {'ReturnValue'}; 
hmmHub.name.modality = hubModalityNames; % just if needed for testing - mainly see below.

% Store in basic structure: -----------------------------------
modelStructure.hmmHub = hmmHub; 
modelStructure.midPFair = midPFair ;
modelStructure.d0 = d0;


%%  ############################################################
%%  ######          Actions :   REPORTING MDPs         #########
%%  ############################################################

%%  V map: allowable policies (of depth T, so T-1 rows). 
%   These are sequences of actions
%   (with an action for each hidden state factor, enumerated in the 
%    last column. T-1 rows bec. actions not taken in last state.)
%  So of the form V(timestep, action-sequence, state-factor) I think
%  REM here stateFNRep = 2;  % No. of state factors: beliefState, reportState                       
%--------------------------------------------------------------------------
V2 = nan(Tsteps2-1,resNRepo,stateFNRep);  % see above.
V2(:,:,1) = 1;  % is actions affecting state factor 1, e.g. beliefHI, and only
                % 'acceptance' - beliefs are changed only due to Hub!
V2(:,:,2) = 1:resNRepo;  %  is actions affecting state factor 2, e.g. reportHI 


%% %% The key dynamics: B map, i.e. Transitions conditioned on action %% %
% controlled transitions: B{u}
%--------------------------------------------------------------------------
% We specify the probabilistic transitions of hidden states for each
% state (as opposed to outcome) factor. Hence actions are themselves
% factorised, state factor F being subect to actions MF :
%  REM from spm_MDP_VB_X: MDP.B{F}(NF,NF,MF)  - transitions among states 
%       under MF control states. So each B{k} refers to *state factor* k,
%       which has Nk states; There must be a grand total of Mk 
%       actions / control states for these, the 'pages' of B{k}, not
%       all of which might be available under V for the current time-step.
%--------------------------------------------------------------------------
% B2{1} : transitions over factor 1, e.g. as-believed-HI : There is only the
%        the action 'accept fate', the pt can't change this 
%        here, only the reporting of it.
B2{1} = zeros(resNRepo,resNRepo,1);
B2{1}(:,:,1)= eye(resNRepo);   
% B{2} is the choice of reported level of e.g. HI
B2{2} = zeros(resNRepo+1,resNRepo+1,resNRepo);
% First, from Initial state which is resolN+1, valid report
% actions lead to the respective report state:
for contrS=1:resNRepo
   noisy = allP.noiseFloor * ones(resNRepo,1); 
   noisy(resNRepo+1) = 0;  % Can't go to these with a valid report action
   noisy(contrS) = 1;    
   noisy = noisy / sum(noisy); 
   B2{2}(:,resNRepo+1,contrS) = noisy;
end
% *from* a valid report control state all actions lead to same
for contrS=1:resNRepo
   % for clarity, first from valid report states:
   B2{2}(1:resNRepo,1:resNRepo,contrS) = eye(resNRepo); 
end


%% %%% outcome probabilities: A template(likelihood map) %%%%%%%%%%% %
%--------------------------------------------------------------------%
% REM from spm_MDP_VB_X: MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes 
%     given hidden states. G must stand for the *outcome factors*, whereas
%     F must stand for the *state factors*, N1 being the number of states
%     in state factor 1 etc. 
% Here specify the probabilistic mapping from hidden states
% to observations (a.k.a. outcomes):
% A2{1} : Informational, (eg HI)-report stage
% A2{2} : How well Intent is to be reported.
% Their form is A{G}(observation, trueState, level-to-report)
% See Emotion_learning_model... A map as example. 
%--------------------------------------------------------------------%

% A2{1} : Informational, (eg HI)-report stage
%             (row)               (col)              (page)   
%      observed-e.g. HIreport  eg as-believed-HI  eg HIreport      
A2{1} = zeros( resNRepo+1,        resNRepo,        resNRepo+1 ); 
% Report outcomes from initial report state, resolN+1, 
% are all 'null', aka noReportMade, only coincidentally resolN+1:
A2{1}(resNRepo+1,:,resNRepo+1) = 1;       % was: A2{1}(resolN+1,:,:,resolN+1,:) = 1;
% Report outcomes from valid report states are as themselves, 
% whatever the true state:
for toRep = 1:resNRepo
    A2{1}(toRep,:,toRep) = 1;
end
%%   A2{2} : How well Intent is to be reported.
%    Done trial-by-trial, see below.   
%            (row)         (col)        (page)                  
%          report quality  attribBelief  attribReport       
%  allP must include, for use here: desBias, resNRepo, Ucor, corLevN, resNRepo :
A2{2} = ARepClassifSerDict(allP);                                                     
% Test with e.g. squeeze(A2{2}(:,1,:)) to see correctness levels for all Hrep when Hbel is 1

%% Reporting MDP priors: (utility) C, aka goal map ------------
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of log
% probabilities over outcomes. 
%    REM from spm_MDP_VB_X: 
%    MDP.C{G}(O,T) - (log) prior preferences for outcomes (modality G),
%    at each of T timesteps or moves.
%    C factors indexed by G must corresp. to A.
% 
%--------------------------------------------------------------------------

% Outcome factor 1 is informational (e.g. HI) report:
C2{1} = zeros(resNRepo+1,Tsteps2 );  % No preferences here!

% Outcome factor for desirability to report what one believes is best: 
%   allP must include, for use here: desCorr, corLevN, Tsteps2
C2{2}  = CRepClassifRepDict(allP); 

%% Element, or unit, of each of the reporting MDPs :  -------

  % We will make working variables mdpHI and mdpSI, which we will populated for every
  % iteration with the correct bits of hmmHub, and then we will use this loop
  % to construct the trial sequence for MDP2 out of the mdpHIs
  %% V, B, C, A, alpha, eta, omega, tau will be the same for all reporting MDPs:
  mdpHI.B = B2;       mdpHI.C = C2;       
  mdpHI.V = V2;    % Consider mdpHI.U = V2;     mdpHI.T = 1;   -- if Ryan OKs ...
  mdpHI.alpha = alphaPrec;
  mdpHI.eta = allP.tiny;         mdpHI.omega = 1-allP.tiny;  % } essentialy flags.
  % Pedantic - the time constant for gradient descent in spm_MDP_VB_XI :
  mdpHI.tau = hmmHub.tau;
  mdpHI.A = A2;  
  
  mdpSI = mdpHI;                      mdpFair = mdpHI;
  
  % States ------------------------------------------------------------------
  mdpHI.s = zeros(stateFNRep,1);      mdpSI.s = mdpHI.s;
  
  % True hidden states of the fairness-reporting MDP. The first row is a 
  % notional true hidden state to report, which shouldn't really play a role.
  mdpFair.s = zeros(stateFNRep,1);           
  mdpFair.s(1,1) = corePBin;         
  mdpFair.s(2,1)   = resNRepo+1;   
  
  % Naming ------------------------------------------------------------------- 
  % All this needed as there's a naming bug in spm_MDP_VB_X*.m
  HIModalityNames = {'whichHIreport','HIreportQual'}; 
  mdpHI.name.modality = HIModalityNames; % see below.
  mdpHI.Aname = {'InfoHIreport','QualityHIreport'};  
  mdpHI.Bname = {'HIbelief','HIreport'};
  SIModalityNames = {'whichSIreport','SIreportQual'}; 
  mdpSI.name.modality = SIModalityNames;
  mdpSI.Aname = {'InfoSIreport','QualitySIreport'};  
  mdpSI.Bname = {'SIbelief','SIreport'};
  fairModalityNames = {'whichFairRep','fairnessRepQual'}; 
  mdpFair.name.modality = fairModalityNames; 
  mdpFair.Aname = {'InfoFairRep','QualityFairRep'};  
  mdpFair.Bname = {'fairnessBelief','fairnessReport'};     


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put together the sequence of trials and solve the model: %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hub: using 'deal' -------------------------------------------------------
[hmmHub(1:trialN)]    = deal(hmmHub);   % create structure array

hmmHub  = spm_MDP_VB_XI(hmmHub);


%% Now do each trial and modality separately! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for trN = 1:trialN
  % Copy from element prototypes:
  mdpH = mdpHI;          mdpS = mdpSI;      mdpF = mdpFair; 
  
  % For the d map of the Fairness reporting MDP, we must use the posterior a 
  % of the PREVIOUS trial, or the very starting one:
  if trN > 1
     mdpF.d = hub2spoke_dFair(hmmHub(trN-1).a{1}, allP);
  else
     mdpF.d = hub2spoke_dFair(hub_a0{1}, allP);   
  end
  
  % true states for the attribution MDPs. REM columns are within-trial timepoints, 
  % rows are state factors. So top row is HItrue, etc.
  mdpH.s(1,1)   = hmmHub(trN).s(1,1);       mdpS.s(1,1)   = hmmHub(trN).s(2,1);
  mdpH.s(2,1)   = resNRepo+1;               mdpS.s(2,1)   = resNRepo+1;
  
  % Now adjust the attribution d maps according to the updated hmmHub. Here, we need 
  %    to know the posterior d and a for the hub MDP, to see what attributes the new
  %    beliefs about a and d amount to. To derive what the 'hi' and 'lo' states
  %    now amount to, we'll  use the portion of aHub that maps 
  %    from hi-lo HI and SI states to returns.
  [mdpH, mdpS] = hub2spoke_dAttr(mdpH, mdpS, hmmHub(trN), allP, pFHS, attrV, iLo, iHi, iHS); 
  % Test with e.g. squeeze(mdpHI.A{2}(:,3,:)) to see correctness levels for all Hreport when Hbelieved is 3

  
  %% ~~~~~~~~~~~~~~~ Solving and Storing the Attibuting MDPs ~~~~~~~~~~~~~~~~~~~~~
  
  % Store the basic structure for output:
  if trN == 1; modelStructure.mdpHI = mdpH; modelStructure.mdpSI = mdpS;  end

  % Now solve and store:
  mdpH  = spm_MDP_VB_XI(mdpH);             mdpS  = spm_MDP_VB_XI(mdpS); 
  mdpF  = spm_MDP_VB_XI(mdpF); 
  
  if trN == 1
      MDPH(1:trialN) = deal(mdpH);  % }  These store first trial and do memory prealloc.
      MDPS(1:trialN) = deal(mdpS);  % }
      MDPF(1:trialN) = deal(mdpF);  % }
   else
      MDPH(trN) = mdpH;     MDPS(trN) = mdpS;    MDPF(trN) = mdpF; %#ok<*SAGROW>
  end
 
  clear mdpH mdpS mdpF; % inefficient but safe ...
  
end

%% Tidy and MDPs for output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%. Tidy up and store the Hub -------------------------------------
% There seems to be a bug naming the modalities within the above, 
% fixed roughly  here:
for k=1:length(hmmHub);  hmmHub(k).label.modality = hubModalityNames; end
MDPs.hub = hmmHub;   % Store 'learning MDP' itself
% Now store d in array form:
MDPs.dHub = {};                    d1 = mdp2arr(hmmHub,'d');
MDPs.dHub{1} = [d0{1} d1{1}];      MDPs.dHub{2} = [d0{2} d1{2}];
% a, the likelihood concentration parameters, in array form:
MDPs.aHub = mdp2arr(hmmHub,'a'); 

% Store 'Reporting MDPs' ------------------------------------------
MDPs.Hattr = MDPH;           MDPs.Sattr = MDPS;       MDPs.fairRep = MDPF; 
% Separately store actions, MDP.u(F,T - 1). 
%  cf. in Step_by_step* tutorial:
% DCM.Y= {MDP.u}; % include the actions made
% Rem: % MDP.u(F,T-1) - vector of actions - for each hidden factor F and for timesteps 1 ... T-1
Resp.Hattr =   {MDPs.Hattr.u};                Resp.Sattr = {MDPs.Sattr.u}; 
Resp.fairRep = {MDPs.fairRep.u};  
% Store observations: 
% Rem: MDP.o(G,T) - matrix of outcomes - for each outcome modality G and timestep T 
% G for the hub are just the fair / unfair returns.
Inp.hub = {MDPs.hub.o}; 
% G for reports are information 'where we are' oucome and 'report quality' outcome.
Inp.Hattr   = {MDPs.Hattr.o};                Inp.Sattr  = {MDPs.Sattr.o}; 
Inp.fairRep = {MDPs.fairRep.o}; 

%% ~~~~~~~~~~~  Plots and demos ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% First a reminder of the learning MDP, ie. hmmHub 
if toPlot > 0
    % Core HMM ------------------------------------------------------------
    spm_figure('GetWin','Hub (learning): First trial'); clf
    % spm_MDP_VB_trial(MDP(1));
    spm_MDP_VB_trial(MDPs.hub(1));

    spm_figure('GetWin','Hub (learning): Last trial'); clf
    spm_MDP_VB_trial(MDPs.hub(trialN));

	% fairness reporting -------------------------------------------------
    spm_figure('GetWin','fairness reporting: First trial'); clf
    % spm_MDP_VB_trial(MDP(1));
    spm_MDP_VB_trial(MDPs.fairRep(1));

    spm_figure('GetWin','fairness reporting: Last trial'); clf
    spm_MDP_VB_trial(MDPs.fairRep(trialN));

    spm_figure('GetWin','fairness reporting: All, trial-by-trial'); clf
    spm_MDP_VB_game(MDPs.fairRep);
      
    disp(' ');
    disp('~~~~~~~~~ Notes to help understand single trial plots: ~~~~~~~~~~~')
    disp('Top left:  Greyscale, strength of post. belief in state for each factor.')
    disp('           Darker is higher prob. Cyan dots are the true underlying state')
    disp('Top right: Belief in which policy will be chosen, darker=surer.')
    disp('           Cyan dot is actually chosen policy.')
    disp('');
    disp('Middle L : Allowable policies, if non-trivial. X axis is the time steps')
    disp('           i.e. the rows of V, y axis the poss. actions, with darkness ')
    disp('           here just the label number of the action - NOT about beliefs etc.')
    disp('Middle R (if present) : evolution of the posterior disribution over policies')
    disp('           over time, from left o right. ')
    disp('');
    disp('Bottom  L: preferences, in the sense of C map priors, over each factor')
    disp('           w. fact.1 uppermost, and darker=more desirable.')
    disp('           [labels appear wrong - factors in serialDictator03 should be')
    disp('             top: info about where we are; ')
    disp('             mid: preferences about reporting accurately, ~[-1,2,3,-1]''')
    disp('                  for wrong,approximate,correct, indifferent/canNotBeBothered')
    disp('             bottom: returns, ~[-2 2 0]'' for unfair,fair,indifferent')
    disp('           Cyan dots are the actual observations/outcomes per timestep; ')
    disp('           in serialDictator01 in the first timestep, before taking action to')
    disp('           report, the pt sees that they are in the undesirable reporting state.')
    disp('Bottom  R: Simulated phasic DA (black) and expected precision (cyan line)')
    disp(' ')
    disp(' ');
    disp('~~~~~~~~~ Notes to help understand muliple trial plots: ~~~~~~~~~~~')
    disp('Topmost  : Circles: true init. state along each factor, from mdp.s, colours are codes:')
    disp('  Modulo numeric <-> {''r.'',''g.'',''b.'',''c.'',''m.'',''k.''}; ')
    disp('  Rows of circles corresp. to factors, first factor TOP-MOST, ')
    disp('  so here trueHI, e.g. 4=cyan; trueSI e.g. green=2; and')
    disp('  initial fairnessLevelReportState 5=magenta')
    disp('2nd down : Observations / outcomes for each factor, from mdp.o ')
    disp('  TOP-MOST is LAST factor here :(, fairness of split, 3=null/not present in timestep, 1=UNFAIR')
    disp('  then how correctly fairness beliefs reported, corrLevN e.g. 3 is the most desirable at timestep 2,')
    disp('  BOTTOM is observation of fairness report made, resolN+1=no report at this timestep')
    disp('           ')
    disp('3rd down : ERP ')
    disp('4th down : DA ')
    disp('5th down : Post. beliefs about states, the y axes bizarelly seem to go HIGHER DOWN! ')
    disp('            ')
    disp('...         ')
end

%--------------------------------------------------------------------------
 if toPlot > 1
    spm_figure('GetWin','Attrib. Harm Intent, trial-by-trial'); clf  % for some reason clf not needed here??
    spm_MDP_VB_game(MDPH);
 end
 if toPlot > 1
    spm_figure('GetWin','Attrib. Selfish Intent, trial-by-trial'); clf  % for some reason clf not needed here??
    spm_MDP_VB_game(MDPS);
 end

%--------------------------------------------------------------------------


if toPlot > 2 % behav. responses for 1st and last trial of Harm Intent report
    spm_figure('GetWin','Attrib. Harm Intent: First trial'); clf
    % spm_MDP_VB_trial(MDP(1));
    spm_MDP_VB_trial(MDPH(1));

    spm_figure('GetWin','Attrib. Harm Intent: Last trial'); clf
    spm_MDP_VB_trial(MDPH(trialN));
end
if toPlot > 2  % behav. responses for 1st and last trial of Selfish Intent report
    spm_figure('GetWin','Attrib. Selfish Intent: First trial'); clf
    spm_MDP_VB_trial(MDPS(1));

    spm_figure('GetWin','Attrib. Selfish Intent: Last trial'); clf
    spm_MDP_VB_trial(MDPS(trialN));
end


if toPlot > 3
    %% Demo of neural dynamics ------------------------------------------------
    % responses to chosen option - 1st trial, incl phase-precession

    %--------------------------------------------------------------------------
	spm_figure('GetWin','Attributing: trial 1 neural responses'); clf
	% spm_MDP_VB_LFP(mdpHI(1),[2 3;3 3],1);
	spm_MDP_VB_LFP(MDPH(1),[1 2;2 2],1);
 
	spm_figure('GetWin','Attributing: last trial neural responses'); clf
	spm_MDP_VB_LFP(MDPH(trialN),[1 2;2 2],1);
   
	% illustrate phase-amplitude (theta-gamma) coupling
    %--------------------------------------------------------------------------
	spm_figure('GetWin','Attributing: all trials neural responses'); clf
	spm_MDP_VB_LFP(MDPH(1:trialN));
    
end  % plotting options

disp('Trial:');
disp(1:trialN);
disp('Returns delivered, 1=unfair, 2=fair:');  
R = mdp2arr(MDPs.hub,'o'); disp(squeeze(R(end,end,:))');

return; 




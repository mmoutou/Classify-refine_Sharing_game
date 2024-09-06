# Classify-refine_Sharing_game

Resources for https://www.biorxiv.org/content/10.1101/2023.05.20.541280v1 and related articles.
To run this code, you will need to clone and make available the functions in 
https://github.com/mmoutou/mmutils and to install SPM from https://www.fil.ion.ucl.ac.uk/spm/software/download/

To see the basic structure of the model, you may wish to run:

I. Demonstration, with plots, of a agent beliefs and learning over a single block
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Use function serialDictator09a.m . To do this, proceed as follows:
1. Copy the present set of files to your computer, under a new folder. You can do this either
2. by downloading and unzipping a zip file, or syncing between this github folder and a ‘local git’ folder.
3. Do the same for the utilities need, mmutils from the same github account, under a different folder in your computer.
4. Add both the above to your matlab path, navigate to 1. above, and open  serialDictator09a.m in your matlab editor.
   Read through comments starting from line 172, and place breakpoints at lines 206, 227, 241, 255, and 284. 
6. Run the following lines at your matlab command line:
close all; 
otherp=[0.5,0.6]; 
selfp= [0.4,0.7, 1,    1/5,    1,    6,  4,0.5,0.1,  0.5,   0]; 
[MDPs,modStruc,Inp,Resp] = serialDictator09(selfp,otherp,1);  

This will first ‘break’ at key points, where you can inspect the elements of the ‘hub’, or learning MDP: 
A, the likelihood map, p(outcome | state) ; 
B, the Transition map p(next state | current state, action ) 
C, Goal priors map ( preferences for outcomes - somewhat analogous to rewards in reinforcement learning ) and 
D - prior probabilities over states. 

5. Press the ‘continue’ button in matlab (double green triangle). This go through the whole demonstration,
6. and produce an ouptut including graphic serialDictator09aDemo.jpg in this github. 
The graphics demonstrate how responses about the probability of fairness are produced
out of the underlying beliefs, accumulated in the ‘hub’. The black horizontal lines in
the top left plot of  serialDictator09aDemo.jpg are the participant’s chosen responses (actions),
while the plot labelled ‘fairness beliefs’ is a grayscale-representation of underlying
belief distributions for each of 12 trials in this block. 

Please refer to the matlab command line output and standard tutorials of active inference 
software for more details about the plots, but researchers who encounter difficulties running 
this code in their own MATLAB systems can contact the corresponding author of the preprint in question, 
https://www.biorxiv.org/content/10.1101/2023.05.20.541280v1.full.pdf+html  ,
or support with this software.

II. Demonstration of full likelihood function 
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

To see how the full likelihood function operates, you may want to step through spm_mdp_L_vi.m (this is for the winning model version). This includes, apart from the basic structure of each block demonstrated above, how learning takes place between blocks and how parameters that relate to ethnicity bias may change the basic parameters.

% First, load a dataset, extract the data for participant 1, block 2 and load the model structure:
load('attrssriaAcIn09.mat'); mapMod.Inp = expD09{1}.partner{2}.Inp; mapMod.Resp = expD09{1}.partner{2}.Resp; load('modStruc09b.mat'); mapMod.mdpStruc = modStruc;

% Then, run the likelihood function:
mL = spm_mdp_L_vi(mapMod.mdpStruc.indexP,mapMod.mdpStruc, mapMod.Inp,mapMod.Resp,1); 

The output structure mL will then contain a detailed record, whose first few items are
      lFair:  Likelihood of 'fair split' data for every trial
      
        lHI:    Likelihood of 'Harm intent' data for every trial
        
        lSI:  Likelihood of 'Self interest' data for every trial
        
        sum:  Sum log likelihood
        
      pDens: sum posterior density over parameters
      
      Lfair: Sum-log-llikelihood for 'fair split' responses
      
      LAttr: Sum log-likelihood for 'attribution' responses.

III. Example of model-fitting through an adaptive grid algorithm 
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

To see an example of model-fitting, run the following function on the matlab command line:

optP09qB1 = attrssri_Grid09q_vi(2,'b','exp_B_Fit09qPt1to66_inclID.mat',12,5,16,0)

Here, the arguments stand for:
2 : This means the second block of blockSize=16 participants, ie starting from pt 17
'b' : This is the follow-up wave of data.
12 : Each grid sweep will have 12*2+1 grid points.
5 :  Measure of fine-ness of grid.
16 :  blockSize=16 participants
0 : Do not run in code testing mode.

The output on the screen will show the grid exploration, listing the participant, then the parameter being explored, the current value of that's being tried out (in transformed space), the task block, and the fit measure. Once each 'sweep' is completed, the current best state of fit is displayed. In the example below, an optimum for pH0Tr=3.5385 is found: 
...
Now pt.=17 for exp_b_Grid09q_17to32 at 9   6  13  31
For transf. pH0Tr=2.6308 Other=1  ,  block sum LL: -24.3317
For transf. pH0Tr=2.6308 Other=2  ,  block sum LL: -53.9456
For transf. pH0Tr=2.6308 Other=3  ,  block sum LL: -47.6423
For transf. pH0Tr=2.6308 Other=4  ,  block sum LL: -42.6943
param. val.: 2.6308  sum post. dens. for all Others: -193.3003
Now pt.=17 for exp_b_Grid09q_17to32 at 9   6  13  31
For transf. pH0Tr=3.5385 Other=1  ,  block sum LL: -22.3972
For transf. pH0Tr=3.5385 Other=2  ,  block sum LL: -53.2924
For transf. pH0Tr=3.5385 Other=3  ,  block sum LL: -46.7334
For transf. pH0Tr=3.5385 Other=4  ,  block sum LL: -42.9037
param. val.: 3.5385  sum post. dens. for all Others: -190.0215
Now pt.=17 for exp_b_Grid09q_17to32 at 9   6  13  31
For transf. pH0Tr=4.4462 Other=1  ,  block sum LL: -22.3972
For transf. pH0Tr=4.4462 Other=2  ,  block sum LL: -53.2924
For transf. pH0Tr=4.4462 Other=3  ,  block sum LL: -47.1506
For transf. pH0Tr=4.4462 Other=4  ,  block sum LL: -42.9794
param. val.: 4.4462  sum post. dens. for all Others: -190.5231
...
Best pH0Tr so far: 3.5385 MAP density:-190.0215
 All best: 3.538539      3.656207      1.692868         0.809             0      1.003613      2.162702      0.666681       0.34004         -0.04      -0.89949             0     -190.0215     -165.3268      350.6536      361.9729
...

Researchers may wish to step through, from innermost loop out:
- Iteration over blocks, line 282
- Iteration over a single parameter grid, line 238
- Iteration over different parameter dimensions, line 143
- Iteration over participants, line 136

Once the funtion has finished running, a full set of fitted parameters, fit measures and descriptives is contained 
in the structure produced ( optP09qB1 above ) and matlab and csv output files are written to disk.

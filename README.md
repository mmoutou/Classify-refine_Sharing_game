# Classify-refine_Sharing_game

Resources for https://www.biorxiv.org/content/10.1101/2023.05.20.541280v1 and related articles.
To run this code, you will need to clone and make available the functions in 
https://github.com/mmoutou/mmutils and to install SPM from https://www.fil.ion.ucl.ac.uk/spm/software/download/

To see the basic structure of the model, you may wish to run:

I. Demonstration, with plots, of a agent beliefs and learning over a single block    
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Use function serialDictator09a.m . To do this, proceed as follows:
1. Copy the present set of files to your computer, under a new folder. You can do this either by downloading and unzipping a zip file, or syncing between this github folder and a ‘local git’ folder.
2. Do the same for the utilities need, mmutils from the same github account, under a different folder in your computer.
3. Add both the above to your matlab path, navigate to 1. above, and open  serialDictator09a.m in your matlab editor. Read through comments starting from line 172, and place breakpoints at lines 206, 227, 241, 255, and 284. 
4. Run the following lines at your matlab command line:
```
close all; 
otherp=[0.5,0.6]; 
selfp= [0.4,0.7, 1,    1/5,    1,    6,  4,0.5,0.1,  0.5,   0]; 
[MDPs,modStruc,Inp,Resp] = serialDictator09(selfp,otherp,1);  
```

This will first ‘break’ at key points, where you can inspect the elements of the ‘hub’, or learning MDP: A, the likelihood map, p(outcome | state) ; B, the Transition map p(next state | current state, action ) C, Goal priors map ( preferences for outcomes - somewhat analogous to rewards in reinforcement learning ) and D - prior probabilities over states. 

5. Press the ‘continue’ button in matlab (double green triangle). This go through the whole demonstration, and produce an ouptut including graphic serialDictator09aDemo.jpg in this github. 
The graphics demonstrate how responses about the probability of fairness are produced out of the underlying beliefs, accumulated in the ‘hub’. The black horizontal lines in the top left plot of  serialDictator09aDemo.jpg are the participant’s chosen responses (actions), while the plot labelled ‘fairness beliefs’ is a grayscale-representation of underlying belief distributions for each of 12 trials in this block. 

Please refer to the matlab command line output and standard tutorials of active inference software for more details about the plots, but researchers who encounter difficulties running this code in their own MATLAB systems can contact the corresponding author of the preprint in question, https://www.biorxiv.org/content/10.1101/2023.05.20.541280v1.full.pdf+html  , for support with this software.

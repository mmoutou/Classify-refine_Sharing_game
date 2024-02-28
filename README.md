# Classify-refine_Sharing_game

Resources for https://www.biorxiv.org/content/10.1101/2023.05.20.541280v1 and related articles.
To run this code, please also clone and make available the functions in 
https://github.com/mmoutou/mmutils

To see the basic structure of the model, you may want to step through serialDictator09a.m to see how the 
basic elements of the MDP (active inference) model are defined 
A - likelihood map, p(outcome | state)
B - Transition map  p(next state | current state, action )
C - Goal priors map ( preferences for outcomes - somewhat analogous to rewards in reinforcement learning ) 
D - prior probabilities over states

To see the full likelihood function, you may want to step through spm_mdp_L_vi.m (and versions other than vi).
This includes, apart from the basic structure, how learning takes place between blocks and how
parameters that relate to ethnicity bias may change the basic parametes.

To see how the models were fitted using an adaptive grid, see matlab functions like 
attrssri_Grid09q.m 
This is the one for the winning model, and is the best commented one.


function [natP, transP, genP, arrP, hd] = genPar_vi_09b(n,fractionalMean,fractionalSD,lo,hi)
% genPar_vi_09b(n,fractionalMean,fractionalSD,lo,hi)
%   n samples are generated from the standard beta distro with 
%   M and SD as per fractionalMean etc., i.e. approx. fraction of the 
%   interval lo - to - high. With defaults for selialDictator06 and logLik fn ...ii 
% For a vaguely broad 'reasonable' distro, [natP,transP,genP,arrP,hd]= genPar_vi_09b(50).
% For a very narrow range, set lo to be the informative one and don't give hi:
%
% REM indexP for betasc: {'pH0','pS0','dInitEv','aInitEv','initEvRat','alphaPrec','wH','wS','w0'}  ;
% modelStructure.priPar=[[1.01 ,1.01,  1.01,     1.01,    1.01,       1.01   ,   10,  10,  10]; ...  % A
%                        [1.01 ,1.01,    2,      2,         2,         2,        10,  10,  10]; ...  % B
%                        [ 0,    0,      0,      0,         0,         0,       -46  -46, -46]; ...  % lo for betasc
%                        [ 1,    1,     100,    100,       100,       100,       46   46   46]];     % hi for betasc
%            %                         <- max        at 1                >       <   SD ~10  >

try fM = fractionalMean; catch, fM=[];  end
try fSD = fractionalSD;  catch, fSD=[]; end
try lo;    catch, lo=[]; end
try hi;    catch, hi=[]; end

hd = {'pH0','pS0','dInitEv','aInitEv','initEvRat','alphaPrec','wH','wS','w0'}

% Now more specific for spm_mdp_L_ii
if  isempty(fM) 
  % see rbetascMS for how 'fractional' below works ...
  % FRACTIONAL mean of  'pH0','pS0','dInitEv','aInitEv','initEvRat','alphaPrec','wH','wS','w0'
  fM =                  [0.2, 0.35,   0.35,    0.65,       0.5,        0.35,   0.5,  0.5, 0.5 ]; 
  % was: FRACTIONAL mean of pH0  pS0   dInit  alphaPrec    wH    wS    desBias  dInitRat
  % fM =                   [0.2, 0.35, 0.35,   0.35,      0.5,  0.5,    0.5,      0.5      ]; 
end
if isempty(fSD)
  % see rbetascMS for how 'fractional' below works ...
  % FRACTIONAL SD of 'pH0','pS0','dInitEv','aInitEv','initEvRat','alphaPrec','wH','wS','w0'
  fSD =             [ 0.1,  0.1,  0.15,       0.15,    0.15,      0.15,     0.15, 0.15, 0.1 ];
  % was: FRACTIONAL SD of   pH0  pS0   dInit  alphaPrec    wH    wS    desBias  dInitRat   
  %    fSD =              [ 0.1, 0.1,  0.15    0.15      0.15   0.15    0.1       0.15 ];
  warning('SDs not given, hence SD, lo and hi all given default values'); 
  % was: lo vs hi of 'pH0','pS0','dInitEv','aInitEv','initEvRat','alphaPrec','wH',  'wS',    'w0'
  lo  =              [0.05,0.05,   2,         2,      0.25,          2,      5.999,  5.999   0.4999 ]; 
  hi  =              [0.95,0.95,   5,         5,      1.75,          15,     6.001,  6.001,  0.5001];  
  % was: lo vs hi of  pH0  pS0   dInit  alphaPrec    wH    wS    desBias  dInitRat   
  % was: lo  =       [0.05,0.05,   2,      2       5.999   5.999   -0.001    0.25 ]; 
  %      hi  =       [0.95,0.95,   5,     15,      6.001,  6.001,   0.001    1.75];  
end
% The following caters for the possibility that we want essentially all params the same:
if isempty(hi)
    warning('hi reset to get narrow intervals')
    hi = lo .* 1.001;
end

arrP = rbetascMS(n,fM,fSD,lo,hi);
natP={};          
transP={};
for ptN = 1:n
    % in native space:
    natP{ptN}.pH0       = arrP(ptN,1);
    natP{ptN}.pS0       = arrP(ptN,2);
    natP{ptN}.dInitEv   = arrP(ptN,3);
    natP{ptN}.aInitEv   = arrP(ptN,4);
    natP{ptN}.initEvRat = arrP(ptN,5);
    natP{ptN}.alphaPrec = arrP(ptN,6);
    natP{ptN}.wH        = arrP(ptN,7);
    natP{ptN}.wS        = arrP(ptN,8);
    natP{ptN}.w0        = arrP(ptN,9);
    
    % transformed space version:
    transP{ptN} = nat2tr_mdp_L_vi(natP{ptN});
end

% Copy of generative macroparameters:
genP.pH0 = [fM(1),fSD(1),lo(1),hi(1)];
genP.pS0 = [fM(2),fSD(2),lo(2),hi(2)];
genP.dInitEv = [fM(3),fSD(3),lo(3),hi(3)];
genP.aInitEv = [fM(4),fSD(4),lo(4),hi(4)];
genP.initEvRat = [fM(5),fSD(5),lo(5),hi(5)];
genP.alphaPrec = [fM(6),fSD(6),lo(6),hi(6)];
genP.wH = [fM(7),fSD(7),lo(7),hi(7)];
genP.wS = [fM(8),fSD(8),lo(8),hi(8)];
genP.w0 = [fM(9),fSD(9),lo(9),hi(9)];

return;


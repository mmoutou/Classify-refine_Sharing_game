function A2 = A2LearnRepDict(allP, midPFair)
% A2LearnRepDict - How well chance-of-fair-return is to be reported.
% midPFair is optional argument, to avoid re-calculating it if already available.
%
%   allP must include, for use here: wH, wS, w0, resolN, corLevN .

wH = allP.wH;           wS = allP.wS;                w0 = allP.w0; 
resolN = allP.resolN;   corLevN = allP.corLevN;

try 
   if isempty(midPFair);  midPFair = midpfair([wH,wS,w0],resolN); end
catch
   midPFair = midpfair([wH,wS,w0],resolN);
end

ubGrid   = (1:resolN)/resolN;

% A2 : How well chance-of-fair-return is to be reported.
%             (row)        (col)   (page)  (issue)                  
%        report quality   trueHI   trueSI  [p(fair)_level,init,final]  
A2 = zeros(corLevN+1,   resolN,  resolN, resolN+2); 

% At first and last report states, all other state components 
% return indifferent/neutral:
A2(corLevN+1,:,:,resolN+1:resolN+2) = 1;         

% Now calc. the resolN 'issues' for time 2 on the basis of the trueHI trueSI, 
% e.g. page A{2}(3, HIlev, SILev, 5 ) would be 1 if state-to-report=5,
% corresponds with HIlev and SIlev. I.e., if HILev and SIlev centres
% combine to give frair return py within lev 5.
%    To have several levels to account for big or small misattributions. Say we
% call diff of 0 OK, up to 1 lev wrong OKish, the rest just wrong.


for kS = 1:resolN   % will fill in A{2}(corLevN-1:corLevN,kH,kS,:)
                    % CHECK WITH e.g.  squeeze(A{2}(1,2,3,:))' to see
                    % which control states, in the last col., correspond
                    % to 'wrong' report (=1) of py of fair offer if HI=2, SI=3.
  for kH = 1:resolN
     % most will be wrong, 
     % so set to wrong etc. and correct later:
     A2(1,kH,kS,1:resolN) = 1;
     % grid cell that contains this is the OK cell: 
     kOK =   find((ubGrid - midPFair(kH,kS)) > 0, 1 );  % ,1 finds the first element in question.
     A2(1      ,kH,kS,kOK) = 0;  % % set Py of 'totally wrong' default to 0 ...
     A2(corLevN,kH,kS,kOK) = 1;  % row corLevN, is the 'correct'
     % One below best, if valid:
     kOKlo = kOK-1;
     if kOKlo > 0
       A2(1,        kH,kS,kOKlo) = 0;  % set Py of 'totally wrong' default to 0
       A2(corLevN-1,kH,kS,kOKlo) = 1;
     end
     % one above best, if valid:
     kOKhi = kOK+1;
     if kOKhi <= resolN  
       A2(1,        kH,kS,kOKhi) = 0;  % set Py of 'totally wrong' default to 0
       A2(corLevN-1,kH,kS,kOKhi) = 1;
      end
  end
end

return;


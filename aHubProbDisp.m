% aHubProbDisp : display how the likelihoods in a are learnt over trials. 

% toPlot = 1;
try
 aHub = MDPs.aHub{1};
catch
    ;
end
aHDim = size(aHub);
T = aHDim(end);
za = sum(aHub,1);
za = squeeze(za);
hiRHubLik = zeros(T,4);
hiRet=2;
for t=1:T
  pa = squeeze( aHub(hiRet,:,:,t)) ./ za(:,:,t); 
  hiRHubLik(t,:) = pa(:);
end

aHubPFig1 = figure;

subplot1 = subplot(1,1,1,'Parent',aHubPFig1);
hold(subplot1,'on');

plot(hiRHubLik);
title('Learning of Lo-Hi HI-SI attribution -> Fair return')
legend('LoLo->Fair','HiLo->Fair','LoHi->Fair','HiHi->Fair')


%  % """""""""""""""""""  Create figure for experimental data """"""""""""""""""""""""""""
%   figure1 = figure;
%   
%   %==========================================================================
%   % First Wheel of Misfortune roulette
%   % subplot 1
%   subplot1 = subplot(2,2,1,'Parent',figure1);
%   view(subplot1,[-76.4145251396648 39.7331778523299]);
%   grid(subplot1,'on');
%   axis(subplot1,'ij');
%   % Set the remaining axes properties
%   set(subplot1,'PlotBoxAspectRatio',[4.11250305473932 1.61803398874989 1],'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
%         'YTick',[1 2 3 4 5]);
%   hold(subplot1,'on');
%   width = 0.25; 
% 
%   bar3(exptbelch(:,1:5,1)',width)   % was ,'g')
%   title('Belief hist. lW 90(70)','FontSize',16)
%   xlabel('report n.','FontSize',12)
%   ylabel('Pshock*10','FontSize',12)
%   zlabel('Prob.','FontSize',12)
%   spm_axis tight, axis auto % square
%   % legend({'Go2Win','Go2AvLoss','NoGo2Win','NoGo2AvLoss'})
% 
%   %==========================================================================
%   % Second Wheel of Misfortune roulette
%   % subplot 2
%   subplot2 = subplot(2,2,2,'Parent',figure1);
%   view(subplot2,[-76.4145251396648 39.7331778523299]);
%   grid(subplot2,'on');
%   axis(subplot2,'ij');
%   % Set the remaining axes properties
%   set(subplot2,'PlotBoxAspectRatio',[4.11250305473932 1.61803398874989 1],'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
%         'YTick',[1 2 3 4 5]);
%   hold(subplot2,'on');
%   width = 0.25; 
% 
%   bar3(exptbelch(:,1:5,2)',width)
%   title('Belief hist. ln 25(30)','FontSize',16)
%   xlabel('report n.','FontSize',12)
%   ylabel('Pshock*10','FontSize',12)
%   zlabel('Prob.','FontSize',12)
%   spm_axis tight, axis auto % square
%   
%   %==========================================================================
%   % Third Wheel of Misfortune roulette
%   % subplot 3




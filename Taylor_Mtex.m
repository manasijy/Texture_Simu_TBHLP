%% This program does texture evolution based on taylor method prescribed by Bunge
% It used duel simplex method to solve taylor ambiguity

clear;
close;
%%
cs = crystalSymmetry('432');
sS = slipSystem.fcc(cs);
ori = txtfile2ori();
q = 0;
epsilon = 0.1*tensor(diag([1 -q -(1-q)]),'name','strain');
    % Other equivalent code to create epsilon tensor 
        % er = tensor([1 -0.5 -0.5],'name','strain')
        % epsilon = 0.3*er.diag;
numIter = 1;
progress(0,numIter);
for sas=1:numIter

  % compute the Taylor factors and the orientation gradients
  [M,~,mori] = calcTaylor_mky(inv(ori) * epsilon ./ numIter, sS.symmetrise,'silent');
  
  % For single orientation
  % [M,b,mori] = calcTaylor(inv(ori)*epsilon,sS.symmetrise);
  
  % rotate the individual orientations
  ori = ori .* mori;
  progress(sas,numIter);
end

plotPDF(ori,Miller({0,0,1},{1,1,1},cs),'contourf')
mtexColorbar

%% Other codes to get more plots from calcTaylor function

% oS = axisAngleSections(cs,cs,'antipodal');
% oS.angles = 10*degree;

% mori = oS.makeGrid;
 
% [M,b,eps] = calcInvTaylor(mori,sS.symmetrise);
%  computing Taylor factor: 100%
% plot(oS,M,'contourf')
% %%
 
 
% [M,b,mori] = calcTaylor(inv(grains.meanOrientation)*epsilon,sS);
% 
% %%
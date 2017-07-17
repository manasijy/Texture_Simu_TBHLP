function [] = plotTexture(o,it_list)

n = length(it_list);
r_c = floor(sqrt(n))+1;
cs = o.CS;
% ss = o.SS;
for i=1:n 
    axesPos = subplot(r_c,r_c,i,'align');
    plotPDF(o(:,i),Miller(1,1,1,cs),'parent',axesPos,'contourf','complete');
    % plotPDF(o(:,i),Miller(1,1,1,cs),'all','parent',axesPos,'MarkerSize',3)%for
    % single orientations one dot looks better
    text(0,0,num2str(it_list(i)));
end


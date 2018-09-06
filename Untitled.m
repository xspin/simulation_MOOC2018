% close all
% clear all
% clc
%open('fig_cov_con/figures/Duration to Sink (s).fig')
%open('fig_cov_con/figures/Nodes Linked to Sink (%).fig')

axesObjs = get(gcf, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
h = findobj(axesObjs, 'type', 'line');
xdata = get(h, 'XData');  %data from low-level grahics objects
ydata = get(h, 'YData');
y=ydata;
ydata{1}

return

i = 1;
n=length(ydata{i});
dy = zeros(1,n);
idx=[11,12,14];
dy(idx) = dy(idx)+0.1;
dy(end)=dy(end)+0.3;
ydata{i}+dy

set(h(i), 'YData', ydata{i}+dy)

for i=2:length(y)
    min(y{1}./y{i})
end

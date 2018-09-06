clear all;close all; clc

% main function for Coverage and Link
parameter=struct('interval_x',[0 2000],...       %(m)
                 'interval_y',[0 2000],...       %(m)
                 'interval_speed',[0 15],...     %(m/s)
                 'velocity', 10,...              %(m/s)
                 'interval_direction',[0, 2*pi],...%(rad)
                 'range_transmission', 130,...      %(m)
                 'range_coverage', 60,...           %(m)
                 'sim_time',200,...          %(s)
                 'sim_time_step', 0.1,...    %(s)
                 'max_hops', 2,...
                 'collision_distance', 5,...    %(m)
                 'sink', 0,...
                 'n_nodes',50);
%-----
% mobility models: CBMM/ABMM/random
% clustering algorithm: MCBCA/MPCA/none
% CB/AB + MC/MP
fh = fopen('log.txt', 'a+');
fprintf(fh, '===============%s=============\n',datestr(datetime('now')));
fprintf(fh, '%s\n', evalc('disp(parameter)'));
fclose(fh);
tstart = tic;
MM = {'', 'RWPMM', 'CVBMM', 'CNBMM', 'ABMM'};
CA = {'', 'MPCA', 'MCBCA'};
MM_style = 'o*';
CA_style = 'gbr';

tic
modellabel = {MM{2:end}, CA{3}};
style = {};
% for j=1:3
%     for i=1:2
%         tp = '+ ';
%         modellabel{end+1} = sprintf('%s%s%s', MM{i},tp(1+isempty(CA{j})), CA{j});
%         style{end+1} = sprintf('%s%s-',MM_style(i), CA_style(j));
%     end
% end
times = 1;
stats = {};
nodes = [30, 50, 70, 100, 150];
nodes = [10, 20];
for n=nodes
    parameter.n_nodes = n;
    fh = fopen('log.txt', 'a+');
    fprintf(fh, 'nodes: %d\n', n);
    fclose(fh);
    stat = [];
    for i=1:length(modellabel)
        fprintf('############# nodes: %d, %d\n',n ,i);
        if i<length(modellabel)
            model.mm = MM{i+1};
            model.ca = CA{1};
        else
            model.mm = MM{1};
            model.ca = CA{3};
        end
        [~, st] = UAV_simulate_stat(parameter, model, times, false);
        stat = [stat; st];
    end
    stats{end+1} = stat;
end
save(sprintf('stats(%s).mat',datestr(now,'DD_HH_MM')), 'stats')
fprintf('Saved statistics data\n')

time = toc(tstart);
[~,~,~,H,M,S] = datevec(time/(24*60^2));
fprintf('\nTotal Elapsed Time: %dh %dm %.0fs\n', H, M, S);
pause

%% plot bars
tt = {'Coverage Rate(%)', 'Duration to Sink(s)', 'Nodes Linked to Sink(%)', 'Link Duration(s)', 'Number of Links', 'Number of CHs', 'CH Duration(s)'};
close all
for t=1:length(tt)
    if any(t==[])
        continue;
    end
    figure
    hold on
    ylabel(tt{t});
    y = stat(:,t);
    lgd = modellabel;
    if t>5
        y = y(3:end);
        lgd = {lgd{3:end}};
    end
    if t==1
        y = y.*100;
    elseif t==3
        y = y./parameter.n_nodes.*100;
    end
    xlim([0, numel(y)+1]);
    for k=1:numel(y)
        %bar(j, y(j));
        plot([k k], [0, max(y(k),max(y)/100)], 'LineWidth',20);
        %text(j-0.3, y(j)+0.1, modellabel{j})
    end
    [~, hobj, ~, ~] = legend(lgd);
    set(hobj,'linewidth',5);
end
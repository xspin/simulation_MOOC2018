clear all;close all; clc

parameter=struct('interval_x',[0 2000],...       %(m)
                 'interval_y',[0 2000],...       %(m)
                 'interval_speed',[1 20],...     %(m/s)
                 'velocity', 10,...              %(m/s)
                 'interval_direction',[0, 2*pi],...%(rad)
                 'range_transmission', 100,...      %(m)
                 'range_coverage', 60,...           %(m)
                 'sim_time',500,...          %(s)
                 'sim_time_step', 0.1,...    %(s)
                 'max_hops', 2,...
                 'collision_distance', 5,...    %(m)
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
MM = {'CBMM', 'ABMM'};
CA = {'', 'MPCA', 'MCBCA'};
MM_style = 'o*';
CA_style = 'gbr';
tic
modellabel = {};
style = {};
for j=1:3
    for i=1:2
        tp = '+ ';
        modellabel{end+1} = sprintf('%s%s%s', MM{i},tp(1+isempty(CA{j})), CA{j});
        style{end+1} = sprintf('%s%s-',MM_style(i), CA_style(j));
    end
end
times = 1;
stats = {};
nodes = [30, 50, 70, 100, 150];
% nodes = [10, 20];
for n=nodes
    parameter.n_nodes = n;
    fh = fopen('log.txt', 'a+');
    fprintf(fh, 'nodes: %d\n', n);
    fclose(fh);
    stat = [];
    for j=1:3
        for i=1:2
            fprintf('############# nodes: %d, %d\n',n ,(j-1)*2+i);
            model.mm = MM{i};
            model.ca = CA{j};
            [~, st] = UAV_simulate_stat(parameter, model, times);
            stat = [stat; st];
        end
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
tt = {'Coverage Rate (%)', 'Duration to Sink (s)', 'Nodes Linked to Sink (%)',...
    'Link Duration (s)', 'Number of Links', 'Number of CHs', 'CH Duration (s)'};
close all


metrics = {};
for t=1:length(tt)     % for metrics
    ydata = [];
    for m=1:length(modellabel)    % for different methods
        tdata = [];
        for n=1:length(stats)
            tdata(end+1) = stats{n}(m, t);
        end
        ydata = [ydata; tdata];
    end
    if t==1
        ydata = ydata * 100;
    elseif t==3
        ydata = ydata ./ nodes.*100;
    end
    metrics{t} = ydata;
end

close all
items = [1, 4, 5, 6, 7];
items = 1:7;
for t=items
    figure
    hold on
    box on
    ylabel(tt{t});
    xlabel('Number of Nodes');
    lgd = modellabel;
    mds = 1:length(lgd);
    if t>5
        lgd = {lgd{3:end}};
        mds = 3:length(mds);
    end
    for i=mds
        plot(nodes, metrics{t}(i,:), style{i});
    end
    legend(lgd);
end















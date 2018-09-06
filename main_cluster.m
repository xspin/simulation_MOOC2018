clear all;close all; clc

% main function for Cluster Performance
parameter=struct('interval_x',[0 1500],...       %(m)
                 'interval_y',[0 1500],...       %(m)
                 'interval_speed',[0 15],...     %(m/s)
                 'velocity', 10,...              %(m/s)
                 'interval_direction',[0, 2*pi],...%(rad)
                 'range_transmission', 200,...      %(m)
                 'range_coverage', 50,...           %(m)
                 'sim_time',500,...          %(s)
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
CA = {'', 'HCCA', 'MPCA', 'MCBCA'};
MM_style = 'o*';
CA_style = 'gbr';

tic
modellabel = {CA{2:3}, 'CBMM'};
style = {'kv-','mo-','b>-','gs-','r*-'};
style = {'bo-', 'gv-', 'r*-'};

%% simulation
times = 5;
stats = {};
nodes = [50, 70, 100, 150];
% nodes = [30, 50];
for n=nodes
    parameter.n_nodes = n;
    fh = fopen('log.txt', 'a+');
    fprintf(fh, 'nodes: %d\n', n);
    fclose(fh);
    stat = [];
    for i=1:length(modellabel)
        fprintf('############# nodes: %d, %d\n',n ,i);
        if i<length(modellabel)
            model.mm = MM{2};
            model.ca = CA{i+1};
        else
            model.mm = MM{1};
            model.ca = CA{4};
        end
        [~, st] = UAV_simulate_stat(parameter, model, times);
        stat = [stat; st];
    end
    stats{end+1} = stat;
    time = toc(tstart);
    [~,~,~,H,M,S] = datevec(time/(24*60^2));
    fprintf(' Current Total Elapsed Time: %dh %dm %.0fs\n', H, M, S);
end
date = datestr(now,'DD_HH_MM');
save(sprintf('fig_cluster/stats(%s).mat',date), 'stats');
fprintf('Saved statistics data\n')

time = toc(tstart);
[~,~,~,H,M,S] = datevec(time/(24*60^2));
fprintf('\nTotal Elapsed Time: %dh %dm %.0fs\n', H, M, S);
% pause

%% plot 
tt = {'Coverage Rate (%)', 'Duration to Sink (s)', 'Nodes Linked to Sink (%)',...
    'Link Duration (s)', 'Average Degree', 'Connectivity',... 
    'Number of CHs', 'CH Duration (s)', 'CM Duration (s)', ... 
    'Number of ONs', 'Change Frequency'};

close all

% date = '15_12_31';
% load(['stats(', date, ').mat'])

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
        ydata = ydata ./ nodes .*100;
    end
    metrics{t} = ydata;
end

save(sprintf('fig_cluster/metrics(%s).mat',date), 'metrics')
fprintf('Saved metrics data\n')

close all
% date = '16_04_51';
% load(sprintf('metrics(%s).mat',date));
folder = ['fig_cluster/fig', date];
mkdir(folder);

items = [1, 4, 5, 6, 7];
items = 7:11;
for t=items
    figure
    hold on
    box on
    ylabel(tt{t});
    xlabel('Number of Nodes');
    lgd = modellabel;
    mds = 1:length(lgd);
    for i=mds
        plot(nodes, metrics{t}(i,:), style{i}, 'LineWidth',2, 'MarkerSize', 10);
    end
    legend(lgd);
    saveas(gcf, [folder '/' tt{t} '.fig']);
end
clear all;close all; clc

% main function for Coverage and Link
parameter=struct('interval_x',[0 2500],...       %(m)
                 'interval_y',[0 2500],...       %(m)
                 'interval_speed',[0 15],...     %(m/s)
                 'velocity', 10,...              %(m/s)
                 'interval_direction',[0, 2*pi],...%(rad)
                 'range_transmission', 120,...      %(m)
                 'range_coverage', 50,...           %(m)
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

MM = {'', 'RWPMM', 'CVBMM', 'CNBMM', 'ABMM'};
CA = {'', 'MPCA', 'MCBCA'};
MM_style = 'o*';
CA_style = 'gbr';

modellabel = {MM{2:end}, 'CBMM'};
style = {'kv-','mo-','b>-','gs-','r*-'};


times = 5;
nodes = [30, 50, 70, 100, 150];
%nodes = [30, 150];


%% plot 
tt = {'Coverage Rate (%)', 'Duration to Sink (s)', 'Nodes Linked to Sink (%)',...
    'Link Duration (s)', 'Average Degree', 'Connectivity',... 
    'Number of CHs', 'CH Duration (s)', 'Number of orphan nodes', 'Reaffiliation Frequency'};


% date = '15_12_31';
% load(['stats(', date, ').mat'])

close all
date = '16_04_51';
load(sprintf('metrics(%s).mat',date));


items = [1, 4, 5, 6, 7];
items = 1:6;
for t=items
    figure
    hold on
    box on
    ylabel(tt{t});
    xlabel('Number of Nodes');
    lgd = modellabel;
    mds = 1:length(lgd);
    if t>6
        lgd = {lgd{3:end}};
        mds = 3:length(mds);
    end
    for i=mds
        plot(nodes, metrics{t}(i,:), style{i}, 'LineWidth',2, 'MarkerSize', 10);
    end
    legend(lgd);
end
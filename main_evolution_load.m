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

MM = {'', 'RWPMM', 'CVBMM', 'CNBMM', 'ABMM'};
CA = {'', 'MPCA', 'MCBCA'};
MM_style = 'o*';
CA_style = 'gbr';

modellabel = {MM{2:end}, 'CBMM'};
style = {'kv-','mo-','b>-','gs-','r*-'};

%% plot
tt = {'Coverage Rate (%)', 'Duration to Sink (s)', 'Nodes Linked to Sink (%)',...
    'Link Duration (s)', 'Average Degree', 'Connectivity',... 
    'Number of CHs', 'CH Duration (s)', 'Number of orphan nodes', 'Reaffiliation Frequency'};

figure
hold on
box on
ylabel(tt{1});
xlabel('Simulation Time (s)');
x = 0:parameter.sim_time_step:parameter.sim_time;
n = 10;
t=1:round(length(x)/n):length(x);
load evol(15_17_02).mat % y

for i=1:length(modellabel)
    %y{i} = sims{i}{1}.area_cover_rate;
    plot(x(t), y{i}(t), style{i}, 'LineWidth',2, 'MarkerSize', 10);
%     plot(x, y, style{i}([1,3]), 'LineWidth',2);
%     plot(x(t), y(t), style{i}([1,2]), 'MarkerSize', 10);
end
legend(modellabel);

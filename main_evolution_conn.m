clear all;close all; clc

% main function for Coverage and Link
parameter=struct('interval_x',[0 2000],...       %(m)
                 'interval_y',[0 2000],...       %(m)
                 'interval_speed',[0 15],...     %(m/s)
                 'velocity', 10,...              %(m/s)
                 'interval_direction',[0, 2*pi],...%(rad)
                 'range_transmission', 130,...      %(m)
                 'range_coverage', 60,...           %(m)
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
% fh = fopen('log.txt', 'a+');
% fprintf(fh, '===============%s=============\n',datestr(datetime('now')));
% fprintf(fh, '%s\n', evalc('disp(parameter)'));
% fclose(fh);
tstart = tic;
MM = {'', 'RWPMM', 'CVBMM', 'CNBMM', 'ABMM'};
CA = {'', 'MPCA', 'MCBCA'};
MM_style = 'o*';
CA_style = 'gbr';

tic
modellabel = {MM{2:end}, 'CBMM'};
% modellabel = {modellabel{4:end}};
style = {'kv-','mo-','b>-','gs-','r*-'};

times = 1;
sims = {};
y = {};
for n=50
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
        [sim, ~] = UAV_simulate_stat(parameter, model, times, false);
        y{i} = sim{1}.connectivity;
    end
end
% save(sprintf('sims(%s).mat',datestr(now,'DD_HH_MM')), 'sims')
% fprintf('Saved sim data\n')

time = toc(tstart);
[~,~,~,H,M,S] = datevec(time/(24*60^2));
fprintf('\nTotal Elapsed Time: %dh %dm %.0fs\n', H, M, S);

%% plot
tt = {'Coverage Rate (%)', 'Duration to Sink (s)', 'Nodes Linked to Sink (%)',...
    'Link Duration (s)', 'Average Degree', 'Connectivity',... 
    'Number of CHs', 'CH Duration (s)', 'Number of orphan nodes', 'Reaffiliation Frequency'};

close all
figure
box on
hold on
ylabel(tt{6});
xlabel('Simulation Time (s)');
x = 0:parameter.sim_time_step:parameter.sim_time;
n = length(y{1});
t=1:round(length(x)/n):length(x);

for i=1:length(modellabel)
    plot(x(t), [1,y{i}], style{i}, 'LineWidth',1, 'MarkerSize', 7);
end
legend(modellabel);
return
save(sprintf('evol_con.mat',datestr(now,'DD_HH_MM')), 'y')
fprintf('Saved evolution data\n')
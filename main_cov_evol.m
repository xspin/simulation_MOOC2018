clear all;close all;

parameter=struct('interval_x',[0 1000],...       %(m)
                 'interval_y',[0 1000],...       %(m)
                 'interval_speed',[1 20],...     %(m/s)
                 'velocity', 10,...              %(m/s)
                 'interval_direction',[0, 2*pi],...%(rad)
                 'range_transmission', 100,...      %(m)
                 'range_coverage', 60,...           %(m)
                 'sim_time',500,...          %(s)
                 'sim_time_step', 0.1,...    %(s)
                 'max_hops', 2,...
                 'collision_distance', 5,...    %(m)
                 'n_nodes',100);
%-----
% mobility models: CBMM/ABMM/random
% clustering algorithm: MCBCA/MPCA/none
% CB/AB + MC/MP
fh = fopen('log.txt', 'a+');
fprintf(fh, '===============%s=============\n',datestr(datetime('now')));
fprintf(fh, '%s\n', evalc('disp(parameter)'));
fclose(fh);

MM = {'CBMM', 'ABMM'};
CA = {'', 'MPCA', 'MCBCA'};
MM_style = {'--','-'};
CA_style = 'gbr';
tic;
modellabel = {};
style = {};
for j=1:3
    for i=1:2
        tp = '+ ';
        modellabel{end+1} = sprintf('%s%s%s', MM{i},tp(1+isempty(CA{j})), CA{j});
        style{end+1} = sprintf('%s%s',MM_style{i}, CA_style(j));
    end
end

times = 1;
sims = {};
for j=1:3
    for i=1:2
        fprintf('#############%d\n',(j-1)*2+i);
        model.mm = MM{i};
        model.ca = CA{j};
        [sims{end+1}, ~] = UAV_simulate_stat(parameter, model, times);
    end
end
toc

%% plot bars
tt = {'Coverage Rate (%)', 'Duration to Sink (s)', 'Nodes Linked to Sink (%)',...
    'Link Duration (s)', 'Number of Links', 'Number of CHs', 'CH Duration (s)'};

close all

figure
hold on
ylabel(tt{1});
xlabel('Simulation Time (s)');
x = 0:parameter.sim_time_step:parameter.sim_time;
for i=1:length(modellabel)
    %y = sims{i}{1}.area_cover_rate;
    plot(x, y{i}, style{i}, 'LineWidth',2);
end
legend(modellabel);



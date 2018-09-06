clear all;close all;

parameter=struct('interval_x',[0 2500],...       %(m)
                 'interval_y',[0 2500],...       %(m)
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
                 'n_nodes',150);
%-----
% mobility models: CBMM/ABMM/random
% clustering algorithm: MCBCA/MPCA/none
% CB/AB + MC/MP

MM = {'', 'RWPMM', 'CVBMM', 'CNBMM', 'ABMM'};
CA = {'', 'HCCA', 'MPCA', 'MCBCA'};
model.mm = MM{1};
model.ca = CA{4};
start=tic;
fprintf('Simulation start\n');             
disp(model)
fprintf('  nodes: %d\n',parameter.n_nodes)
sim = MobilitySimulation(parameter, model);
sim = sim.initialize(); 
sim = sim.simulate();
fprintf(' Statistics')
a = tic;
st = sim.statistics()
toc(a)
toc(start)


sim.animate();


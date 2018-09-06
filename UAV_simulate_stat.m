function [sims, stat] = UAV_simulate_stat(parameter, model, times, flag)
    % simulate the model for multipy times
    if nargin<4
        flag = true;
    end
    fh = fopen('log.txt', 'a+');
    fprintf('Start simulatation\n');
    fprintf('  Mobility Model: %s\n', model.mm);
    fprintf('  Clustering Alg.: %s\n', model.ca);
    fprintf('  Repeat for %d times\n', times);
    sims = {};
    stat = 0;
    tstart = tic;
    for t=1:times
        fprintf(' Simulation: %d\n',t);             
        sim = MobilitySimulation(parameter, model);
        sim = sim.initialize(); 
        sim = sim.simulate();
        fprintf('\n stat: ')
        [st, sim] = sim.statistics();
        fprintf('\n')
        stat = stat + struct2array(st);
        time = toc(tstart);
        [~,~,~,H,M,S] = datevec(time/(24*60^2));
        fprintf('Elapsed time: %dh %dm %.0fs\n', H, M, S);
        if t<times
            [~,~,~,H,M,S] = datevec((time*times/t-time)/(24*60^2));
            fprintf('Expected remaining time: %dh %dm %.0fs\n', H, M, S);
        end
        sims{end+1}=sim;
        if flag
            fprintf(fh, '%s|%s|%d: %.2f', model.mm,model.ca, t);
            fprintf(fh, '%.2f,', struct2array(st));
            fprintf(fh, '\n');
        end
        fprintf('%.2f,', struct2array(st));
        fprintf('\n')
    end
    stat = stat/times;
    if flag
        fprintf(fh, '%s|%s| : %.2f', model.mm, model.ca);
        fprintf(fh, '%.2f,', stat);
        fprintf(fh, '\n');
    end
    fclose(fh);
end
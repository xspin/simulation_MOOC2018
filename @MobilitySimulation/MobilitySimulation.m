classdef MobilitySimulation
    %MOBILITYSIMULATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model;
        interval_x;         % [a, b]
        interval_y;         % [a, b]
        interval_speed;     % [a, b]
        velocity;           % v
        interval_direction; % [a, b]
        range_transmission; % D
        range_coverage;     % r
        sim_time;           % t
        sim_time_step;      % t
        max_hops;           % d
        n_nodes;            % N
        sink;
        %fid;                % first UAV id
        %min_cm;
        fis;
        
        nodes;              % [n1, n2, ...]
        trace;              % {t1, t2, ...}
        chs;                % [ni, nj, ...]
        area_coverage;      % area(i, j) = visited_time;
        area_cover_rate;    % [t...]
        area_discrecity;
        
        % metrics
        ch_count;
        connectivity;       % [t...]

    end
    
    methods
        function obj = MobilitySimulation(input, model)
            obj.interval_x = input.interval_x;
            obj.interval_y = input.interval_y;
            obj.interval_speed = input.interval_speed;
            obj.velocity = input.velocity;
            obj.interval_direction = input.interval_direction;
            obj.range_transmission = input.range_transmission;
            obj.range_coverage = input.range_coverage;
            obj.sim_time = input.sim_time;
            obj.sim_time_step = input.sim_time_step;
            obj.max_hops = input.max_hops;
            obj.n_nodes = input.n_nodes;
            obj.sink = input.sink;
            %obj.min_cm = input.min_cm;
            obj.model = model;
            obj.fis = fis_ABMM();
            
        end
        
        % function obj = initialize(obj)
        % function obj = simulate(obj)
        % function obj = statistics(obj)
 
        function obj = record(obj)
            
            % nodes's states
            obj.trace{end+1} = obj.nodes;
            
            % coverage metrics
            for i=1:obj.n_nodes
                if i==1
                    continue
                end
                idx = ceil(obj.nodes(i).p/obj.area_discrecity);
                obj.area_coverage(idx(2), idx(1)) = obj.area_coverage(idx(2), idx(1))+obj.sim_time_step;
            end
            data = reshape(obj.area_coverage, 1, numel(obj.area_coverage));
            obj.area_cover_rate(end+1) = sum(data>0)/numel(data);
            
        end
        
        function plot_coverage(obj)
            data = reshape(obj.area_coverage, 1, numel(obj.area_coverage));
            cov_min = min(data);
            cov_max = max(data);
            cov_avg = mean(data);
            cov_med = median(data);
            cov_var = var(data);
            cov_zero = sum(data==0);
            fprintf('coverage:\n min=%.2f, max=%.2f, med=%.2f, avg=%.2f, var=%.2f\n',...
                cov_min, cov_max, cov_med, cov_avg, cov_var);
            fprintf(' zeros:%d\n', cov_zero);
            
            figure;
            hist(data, 100);
            
            figure;
%             box on;
%             grid on;
            axis([obj.interval_x, obj.interval_y]);
            title('Coverage');
            [x, y] = meshgrid(1:size(obj.area_coverage, 2), 1:size(obj.area_coverage, 1));
            x = x.*obj.area_discrecity - obj.area_discrecity/2;
            y = y.*obj.area_discrecity - obj.area_discrecity/2;
            surf(x, y, obj.area_coverage);
            shading interp
            colorbar
            figure;
            title('Coverage');
            pcolor(x,y, obj.area_coverage);
%             shading interp
            colorbar
            
        end

        function hd = plot(obj, nodes)
            if nargin<2
                nodes = obj.nodes;
            end
            hd = {};
            hold on
            for node = nodes
                color = node.type;
                if color == 'b'
                   color = 'k';
                elseif color=='g'
                    color = [0.7,0.7,0.7];
                end
                % plot nodes
                hd{end+1,1} = plot(node.p(1), node.p(2), 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',color);
                % plot cluster links
                hd{end,2} = plot([0, 0],[0, 0],'-.', 'Color', [0.5, 0.5, 0.5]);
                % plot sink links
                hd{end,3} = plot([0, 0],[0, 0],'b-.');
                % display cluster links' length
                hd{end,4} = text(0,0,'');
                % display cluster heads' ID
                hd{end,5} = text(0,0,'', 'Color', 'r');
                % plot links between clusters to sink
                hd{end,6} = plot([0, 0],[0, 0],'r-.');
            end
            if obj.sink>0
                theta = 0:0.1:2*pi;
                xy = obj.range_transmission*[cos(theta); sin(theta)]+obj.nodes(1).p;
                plot(xy(1,:), xy(2,:), 'r--');
            end
            hold off
        end
             
        function animate(obj, Time)
            flag_sink = false;
            flag_dist = false;
            if nargin<2
                Time = 1:length(obj.trace);
            else
                Time = round(Time/obj.sim_time_step);
                if length(Time)>2
                    Time = Time(1):Time(2);
                end
            end
            figure;
            box on;
            grid on;
            axis([obj.interval_x, obj.interval_y]);
            %title('UAV Mobility Simulation');
            xlabel('X (m)');
            ylabel('Y (m)');
            ht = text(obj.interval_x(1), obj.interval_y(2)-13, sprintf('Time : 0s'));
            hd = obj.plot(obj.trace{1});
            for t=Time
                nodes = obj.trace{t};
                set(ht,'String',sprintf('Time : %.0fs', t*obj.sim_time_step));
                for i=1:obj.n_nodes
                    set(hd{i, 1}, 'XData', nodes(i).p(1), 'YData', nodes(i).p(2));
                    color = nodes(i).type;
                    if color == 'b'
                       color = 'k';
                    elseif color=='g'
                        color = [1,1,1]*0.7;
                        %color = 'b';
                    end
                    if i==1
                        color = 'r';
                    end
                    set(hd{i,1}, 'MarkerFaceColor',color);
                    if nodes(i).type=='g'
                        set(hd{i,2}, 'XData', [nodes(i).p(1), nodes(nodes(i).pid).p(1)], 'YData', [nodes(i).p(2), nodes(nodes(i).pid).p(2)]);
                        if flag_dist
                            set(hd{i,4}, 'Position',(nodes(nodes(i).pid).p+nodes(i).p)/2,'String', sprintf('%.0f', norm((nodes(nodes(i).pid).p-nodes(i).p))));
                        end
                    else
                        set(hd{i,2}, 'XData', [], 'YData', []);
                        set(hd{i,4}, 'String', '');
                    end
                    if nodes(i).type=='b' && false
                        set(hd{i,5}, 'Position',nodes(i).p+[2;2], 'String', num2str(i));
                        
                    else
                        set(hd{i,5}, 'String', '');
                    end
                    if nodes(i).cpid>0
                        set(hd{i,6}, 'XData', [nodes(i).p(1), nodes(nodes(i).cpid).p(1)], 'YData', [nodes(i).p(2), nodes(nodes(i).cpid).p(2)]);
                    else
                        set(hd{i,6}, 'XData', [], 'YData', []);                        
                    end
                    if flag_sink
                        if nodes(i).pid2sink>0
                            set(hd{i, 3}, 'XData', [nodes(i).p(1), nodes(nodes(i).pid2sink).p(1)], 'YData', [nodes(i).p(2), nodes(nodes(i).pid2sink).p(2)]);
                        else
                            set(hd{i,3}, 'XData', [], 'YData', []);
                        end
                    end
                end
                drawnow;
            end
        end
    end
    
end


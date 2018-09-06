classdef Node
    properties
        id;     % 1~N
        type;   % 'w'/'g'/'b'
        p;      % [x; y]
        v = 0;      % velocity
        yaw = 0;    % 
        chid = 0;   % CH ID
        pid = 0;    % pre ID to CH
        hop = 0;    % hop to CH
        pid2sink = 0;
        hop2sink = Inf;
        weight = 0;     % followship weight
        weight_mpca = 0;
        neighbor = [];
        neighbor_cluster = [];
        cluster_member = [];
        cpid = 0;
        force = [];
        pd = 100;
    end
    
    methods
        function node = Node(id, p, type)
            if nargin<1
                return;
            end
            node.id = id;
            node.p = p;
            node.type = type;
        end
        
%         function display(obj)
%             fprintf('\tID:%d, Type:%s, CH ID:%d\n\tPos:(%.2f,%.2f), Yaw:%.2f, Velocity:%.2f\n',...
%                 obj.id, obj.type, obj.hid, obj.p, obj.yaw, obj.v);
%         end
        function y = inbound(obj, bound_x, bound_y, tmp_p)
            if nargin<4
                tmp_p = obj.p;
            end
            y = false;
            if bound_x(1)<tmp_p(1) && tmp_p(1)<bound_x(2) ...
               && bound_y(1)<tmp_p(2) && tmp_p(2)<bound_y(2)
               y = true;
            end
        end
        
        function obj = move(obj, t, bound_x, bound_y)
            %while t>1e-4
            for k=1:5
                tmp_p = obj.p + t * obj.v * [cos(obj.yaw); sin(obj.yaw)];
                if obj.inbound(bound_x, bound_y, tmp_p)
                    obj.p = tmp_p;
                    break;
                else
                    obj.yaw = mod(obj.yaw, 2*pi);
                    if k>=2
                        % obj.yaw = rand()*pi*2;
                        error('fuck')
                    end
                    target_p = [bound_x(2); bound_y(2)]/2;
                    tmp = target_p-obj.p;
                    obj.yaw = atan(tmp(2)/tmp(1))+(tmp(1)<0)*pi;
                    obj.yaw = mod(obj.yaw, 2*pi);
                end
            end
        end
        
        function obj = bewhite(obj)
            obj.type = 'w';
            obj.chid = 0;
            obj.pid = 0;
            obj.hop = 0;
            obj.cluster_member = [];
        end
    end
end
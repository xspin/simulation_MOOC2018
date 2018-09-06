function obj = initialize(obj)
    obj = node_generate(obj);
    obj.chs = [];
    for i=1:length(obj.nodes)
        obj.nodes(i).v = obj.velocity;
        obj.nodes(i).yaw = rand()*pi*2;
    end
    if obj.sink>0
        % node 1 is the sink
        obj.nodes(1).v = 0;
        obj.nodes(1).p = [sum(obj.interval_x); sum(obj.interval_y)]/2;
    end
    % discrete area
    obj.area_discrecity = obj.range_coverage;
    rows = ceil(obj.interval_y(2)/obj.area_discrecity);
    cols = ceil(obj.interval_x(2)/obj.area_discrecity);
    obj.area_coverage = zeros(rows, cols);
    
end

function obj = node_generate(obj)
    d = obj.range_transmission/2+10; %(m)
    p0 = [sum(obj.interval_x); sum(obj.interval_y)]/2;
    n_row = floor(sqrt(obj.n_nodes));
    obj.nodes = Node();
    
    for id=1:obj.n_nodes
        % random position
        w = sqrt(obj.n_nodes)*20;
        %w = obj.interval_x(2)/2;
%         if strcmp(obj.model.ca, 'MCBCA')
%             p0 = [w;w];
%         end
        obj.nodes(id) = Node(id, p0+[-w+rand*w*2; -w+rand*w*2], 'w');
    end
%     if strcmp(obj.model.ca, 'MCBCA')
%         obj.nodes(1).p = [10; 10];
%         obj.nodes(1).yaw = pi/4;
%     end
    return
    
    for id=1:obj.n_nodes
        row = floor((id-1)/n_row)+1;
        col = mod(id-1, n_row)+1;
        obj.nodes(id) = Node(id, p0+[col*d; row*d], 'w');
    end
end
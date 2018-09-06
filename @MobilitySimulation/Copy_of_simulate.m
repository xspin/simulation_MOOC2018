function obj = simulate(obj)
% simulation process
    for time=0:obj.sim_time_step:obj.sim_time
        if mod(time, obj.sim_time/10)==0
            fprintf('%3.0f%%', time/obj.sim_time*100);
        end
        
        % move to next step
        for i=obj.sink+1:obj.n_nodes
            obj.nodes(i) = obj.nodes(i).move(obj.sim_time_step, obj.interval_x, obj.interval_y);
        end
        
        % clear data
        for i=1:obj.n_nodes
            obj.nodes(i).neighbor = [];
            obj.nodes(i).neighbor_cluster = [];
            obj.nodes(i).cluster_member = [];
            obj.nodes(i).cpid = 0;
            obj.nodes(i).pid = 0;
            obj.nodes(i).hop = 0;
            obj.nodes(i).pid2sink = 0;
            obj.nodes(i).hop2sink = Inf;
            if obj.nodes(i).type == 'g'
                obj.nodes(i).type = 'w';
            end
        end
        
        % update neighbors
        for i=1:obj.n_nodes
            for j=i+1:obj.n_nodes
                if norm(obj.nodes(i).p-obj.nodes(j).p) <= obj.range_transmission
                    obj.nodes(i).neighbor(end+1) = j;
                    obj.nodes(j).neighbor(end+1) = i;
                end
            end
        end
        
        % cluster heads broadcast HELLO to its cluster members
        obj = broadcast_hello(obj);
        
        % choose a clustering algorithm
        if strcmp(obj.model.ca, 'MPCA')
            obj = clustering_MPCA(obj);
        elseif strcmp(obj.model.ca, 'MCBCA')
            obj = clustering_MCBCA(obj, time);
        elseif strcmp(obj.model.ca, 'HCCA')
        end
        
        % update the members linked to sink by BFS
        nb_sink = 5;
        queue = {[0, 1]};
        while ~isempty(queue)
            for i=queue{1}(2:end)
                hop = queue{1}(1)+1;
                tmp = [];
                neighbor = obj.nodes(i).neighbor;
                if i==1 && length(neighbor)>nb_sink % limite the neighors of sink
                    dist = sum(([obj.nodes(neighbor).p]-obj.nodes(i).p).^2);
                    [~, I] = sort(dist);
                    neighbor = neighbor(I(1:nb_sink));
                end
                for j=neighbor
                    if j==1 || obj.nodes(j).pid2sink>0
                        continue;
                    end
                    obj.nodes(j).pid2sink = i;
                    obj.nodes(j).hop2sink = hop;
                    tmp(end+1) = j;
                end
                queue{end+1} = [hop, tmp];
            end
            queue(1) = [];
        end
        
        % broadcast between clusters
        obj = cluster_hello(obj);
        
        % update links ruoting to sink
        obj = update_route(obj);
        
        % update physical state of nodes / mobility model
        if strcmp(obj.model.mm, 'CVBMM')
            obj = update_state_CVBMM(obj);
        elseif strcmp(obj.model.mm, 'CNBMM')
            obj = update_state_CNBMM(obj);
        elseif strcmp(obj.model.mm, 'ABMM')
            obj = update_state_ABMM(obj);
        elseif strcmp(obj.model.mm, 'RWPMM')
            obj = update_state_RWPMM(obj);
        end
        if strcmp(obj.model.ca, 'MCBCA')
            obj = update_state_MCBCA(obj);

        end
        pt = 1;
        wt = 10;
        pause_time = pt/obj.sim_time_step; 
        walk_time = wt/obj.sim_time_step;  
        p1 = 1/(pause_time+1);
        p2 = 1/(walk_time+1);
        for i=2:obj.n_nodes
            break
            if obj.nodes(i).v == 0 && rand<p1 % pause
                obj.nodes(i).v = obj.velocity;
                obj.nodes(i).yaw = rand * 2*pi;
            elseif obj.nodes(i).v>0 && rand<p2 % walk
                obj.nodes(i).v = 0;
            end
        end
        
        obj = update_state(obj);
 
        % record current status / metrics
        obj = obj.record();
    end
end

%% 
function myprintf(varargin)
    debug = false;
    if debug
        fprintf(varargin{:});
    end
end

function obj = cluster_hello(obj)
    for ch_id=obj.chs
        for i=[ch_id, obj.nodes(ch_id).cluster_member]
            nb = [i, obj.nodes(i).neighbor];
            nb = nb([obj.nodes(nb).chid]>ch_id);
            for j=nb
                jchid = obj.nodes(j).chid;
                if ~ismember(jchid, obj.nodes(ch_id).neighbor_cluster)
                    obj.nodes(ch_id).neighbor_cluster(end+1) = jchid;
                    obj.nodes(jchid).neighbor_cluster(end+1) = ch_id;
                end
            end
        end
    end
end

function obj = update_route(obj)
    queue = {[0, 1]};
    while ~isempty(queue)
        for ch_id=queue{1}(2:end)
            hop = queue{1}(1)+1;
            tmp = [];
            if ch_id==1
                nb = obj.nodes(ch_id).neighbor;
                neighbor = unique([obj.nodes(nb).chid]);
                if isempty(neighbor)
                    neighbor = [];
                end
            else
                neighbor = obj.nodes(ch_id).neighbor_cluster;
            end
            for j=neighbor
                if j<=1 || obj.nodes(j).cpid>0
                    continue;
                end
                obj.nodes(j).cpid = ch_id;
                obj.nodes(j).hop2sink = hop;
                tmp(end+1) = j;
            end
            queue{end+1} = [hop, tmp];
        end
        queue(1) = [];
    end
end

function obj = broadcast_hello(obj, chs)
    % broadcast Hello from CH node id to its CMs for keep-alive
    % breadth first search
    if nargin<2
        chs = obj.chs;
    end
    for ch_id = chs
        queue = {[0, ch_id]};
        while ~isempty(queue)
            for i=queue{1}(2:end)
                hop = queue{1}(1)+1;
                tmp = [];
                for j=obj.nodes(i).neighbor
                    if obj.nodes(j).type == 'b'||obj.nodes(j).type == 'g'
                        continue;
                    end
                    if obj.nodes(j).chid == ch_id
                        % node j is still in the cluster
                        obj.nodes(ch_id).cluster_member(end+1) = j;
                        obj.nodes(j).type = 'g';
                        obj.nodes(j).pid = i;
                        obj.nodes(j).hop = hop;
                        tmp(end+1) = j;
                    end
                end
                if hop < obj.max_hops
                    queue{end+1} = [hop, tmp];
                end
            end
            queue(1) = [];
        end
    end
end

function obj = update_state(obj)

    % adjust yaw toward less covered area
    a = 0.9;
    for i=1:obj.n_nodes
        break
        target_p = next_area(obj, obj.nodes(i).p);
        tmp = target_p-obj.nodes(i).p;
        yaw = atan(tmp(2)/tmp(1))+(tmp(1)<0)*pi;
        yaw = mod(yaw, 2*pi);
        if abs(obj.nodes(i).yaw - yaw)>pi
            if obj.nodes(i).yaw>yaw
                obj.nodes(i).yaw = obj.nodes(i).yaw-2*pi;
            else
                yaw = yaw-2*pi;
            end
        end
        obj.nodes(i).yaw = a*obj.nodes(i).yaw + (1-a)*yaw;
        obj.nodes(i).yaw = mod(obj.nodes(i).yaw, 2*pi);
        if isnan(obj.nodes(i).yaw)
            error('NAN yaw of Node-%d in toward an area\n',i);
        end
    end
    
    % limite speed
    for i=obj.sink+1:obj.n_nodes
        obj.nodes(i).v = max(obj.nodes(i).v, obj.interval_speed(1));
        obj.nodes(i).v = min(obj.nodes(i).v, obj.interval_speed(2));
    end
    if obj.sink>0
        obj.nodes(1).v = 0;
        obj.nodes(1).p = [sum(obj.interval_x); sum(obj.interval_y)]/2;
    end
end

%% MCBCA
function obj = clustering_MCBCA(obj, time)
% Mobility Control Based Clustering Algorithm
    % white nodes join a nearest cluster
    obj = join_nearest_cluster(obj);

    % white nodes form clusters
    obj = clustering(obj);

    % adjust/change cluster heads
    if mod(time, 10/obj.sim_time_step)==0
        obj = adjust_ch(obj);
    end
end

function obj = clustering(obj)
    min_cm = 2;
    ch = [];
    for i=1:length(obj.chs)
        ch_id=obj.chs(i);
        if length(obj.nodes(ch_id).cluster_member)<min_cm
            for j=obj.nodes(ch_id).cluster_member
                obj.nodes(j) = obj.nodes(j).bewhite();
            end
            obj.nodes(ch_id) = obj.nodes(ch_id).bewhite();
            ch(end+1) = i;
        end
    end
    if ~isempty(ch)
        myprintf(' Delete CH: '); myprintf('%d ', obj.chs(ch)); myprintf('\n');
        obj.chs(ch) = [];
        myprintf(' CHs: '); myprintf('%d ', obj.chs); myprintf('\n');
    end
    id = 0;
    wnb = [];
    for i=1:obj.n_nodes
        if obj.nodes(i).type == 'w'
            nb = obj.nodes(i).neighbor;
            twnb = nb([obj.nodes(nb).type]=='w');
            if (length(twnb)>length(wnb))
                id = i;
                wnb = twnb;
            end
        end
    end
    if length(wnb)<min_cm
        return;
    end
    obj.chs(end+1) = id;
    myprintf(' Add CH: %d\n', id);
    obj.nodes(id).type = 'b';
    obj.nodes(id).cluster_member = wnb;
    obj.nodes(id).chid = id;
    for i=wnb
        obj.nodes(i).type = 'g';
        obj.nodes(i).hop = 1;
        obj.nodes(i).chid = id;
        obj.nodes(i).pid = id;
    end
    myprintf(' CHs: '); myprintf('%d ', obj.chs); myprintf('\n');
end

function target = next_area(obj, p)
    idx = ceil(p/obj.area_discrecity);
    r = idx(2);
    c = idx(1);
    cov = Inf;
    d = 4;
    for rr=max(r-d,1):min(r+d,size(obj.area_coverage,1))
        for cc=max(c-d,1):min(c+d,size(obj.area_coverage,2))
            if obj.area_coverage(rr,cc)< cov
                idx = [cc; rr];
                cov = obj.area_coverage(rr,cc);
            end
        end
    end
    target = idx*obj.area_discrecity;
end
        
function obj = adjust_ch(obj)
    chs = [];
    for i=1:length(obj.chs)
        ch_id = obj.chs(i);
        vec = [cos(obj.nodes(ch_id).yaw-pi/2), sin(obj.nodes(ch_id).yaw-pi/2)];
        dist = vec * ([obj.nodes([obj.nodes(ch_id).cluster_member, ch_id]).p]);
        [~, I] = sort(dist);
        mid = I(floor(length(I)/2)+1);
        if mod(length(dist), 2)==0
            if I(floor(length(I)/2))==length(I) || mid==length(I)
                continue;
            end
        else
            if mid==length(I)
                continue;
            end
        end
        j = obj.nodes(ch_id).cluster_member(mid);
        obj.nodes(j).type = 'b';
        obj.nodes(j).chid = j;
        obj.nodes(j).cluster_member = obj.nodes(ch_id).cluster_member;
        obj.nodes(j).cluster_member(mid) = ch_id;
        obj.nodes(ch_id).type = 'w';
        obj.nodes(ch_id).chid = j;
        obj.nodes(ch_id).cluster_member = [];
        obj.chs(i) = j;
        chs(end+1) = j;
        myprintf(' Exchange CHs: %d <-> %d\n', ch_id, j);
        myprintf(' CHs: '); myprintf('%d ', obj.chs); myprintf('\n');
    end
    obj = broadcast_hello(obj, chs);
end

function obj = join_nearest_cluster(obj)
    for i=1:obj.n_nodes
        if obj.nodes(i).type ~= 'w'
            continue
        end
        tmp = 0;
        for j=obj.nodes(i).neighbor
            if obj.nodes(j).type == 'w'
                continue;
            end
            if tmp==0 || obj.nodes(j).hop<obj.nodes(tmp).hop ...
                || (obj.nodes(j).hop==obj.nodes(tmp).hop &&... 
                    norm(obj.nodes(i).p-obj.nodes(j).p)<norm(obj.nodes(i).p-obj.nodes(tmp).p))
                tmp = j;
            end
        end
        if tmp>0 && obj.nodes(tmp).hop<obj.max_hops
            ch_id = obj.nodes(tmp).chid;
            obj.nodes(i).type = 'g';
            obj.nodes(i).chid = ch_id;
            obj.nodes(i).hop = obj.nodes(tmp).hop+1;
            obj.nodes(i).pid = tmp;
            obj.nodes(ch_id).cluster_member(end+1) = i;
        end
    end
end

function f = virtual_force(obj, type, i, j)
    a;
    kn;
    kv;
    kl;
    Rt = obj.range_transmission;
    Rc = obj.range_coverage;
    %ch = obj.nodes(i).chid;
    yaw = obj.nodes(j).yaw;
    pij = obj.nodes(j).p-obj.nodes(i).p;
    dij = norm(pij);
    eij = pij/dij;
    vij = [sin(yaw), -cos(yaw)];
    hij = dot(pij, vij);
    if hij>0
        vij = -vij;
    else
        hij = -hij;
    end
    if strcmp(type, 'con')
        f = kn*(dij-a*Rt) * eij;
    elseif strcmp(type, 'cov')
        f = kv*(2*Rc-min(hij, 2*Rc)) * vij;
    elseif strcmp(type, 'loo')
        f = kl*eij/dij; 
    end
end

function obj = update_state_CBMM(obj)
    force = zeros(1, obj.n_nodes);
    for i=2:obj.n_nodes
        if obj.nodes(i).type ~= 'w'
            % mobility for ONs
            for j=obj.nodes(i).neighbor
                force(i) = force(i) + virtual_force(obj, 'con', i, j);
            end
        elseif obj.nodes(i).type == 'g'
            % mobility for CMs
            j=obj.nodes(i).pid;
            force(i) = force(i) + virtual_force(obj, 'cov', i, j);
        elseif obj.nodes(i).type == 'b'
            % mobility for CHs
            force(i) = force(i) + virtual_force(obj, 'loo', i, j);
        end
    end
end

function obj = update_state_MCBCA(obj)
    da = 0.6*obj.range_transmission;
    db = 0.8*obj.range_transmission;
    dc = (da+db)/2;
    if 0
        % Update cluster's direction for keeping connection
        k = 1;
        % update white nodes' direction
        for i=obj.sink+1:obj.n_nodes
            if obj.nodes(i).type ~= 'w'
                continue
            end
            p = obj.sim_time_step/0.5;
            if rand>p
                continue
            end
            nb = obj.nodes(i).neighbor;
            if isempty(nb)
                continue
            end
            pij = [obj.nodes(nb).p]-obj.nodes(i).p;
            dij = sqrt(sum(pij .* pij));
            pim = pij .* (1-dc./dij);
            Fij = k*pim .* or(dij>db, dij<da);
            Fi = sum(Fij,2);
            t = rand*2*pi;
            yaw = vang(Fi + [cos(t); sin(t)]);
            if ~isnan(yaw)
                obj.nodes(i).yaw = yaw;
            end
        end
        % update cluster direction based on forces
        for ch_id=obj.chs
            F = [0; 0];
            p = obj.sim_time_step/1;
            if rand>p
                continue
            end
            for i=[ch_id, obj.nodes(ch_id).cluster_member]
                nb = obj.nodes(i).neighbor;
                ch = [obj.nodes(nb).chid];
                nb = nb(ch ~= ch_id);
                if isempty(nb)
                    continue
                end
                pij = [obj.nodes(nb).p]-obj.nodes(i).p;
                dij = sqrt(sum(pij .* pij));
                pim = pij .* (1-dc./dij);
                Fij = k*pim .* or(dij>db, dij<da);
                t = sum(Fij,2);
                if norm(t)>norm(F)
                    F = t;
                end
            end
            t = rand*2*pi;%+ [cos(t); sin(t)]
            F = F+[cos(obj.nodes(ch_id).yaw); sin(obj.nodes(ch_id).yaw)];
            yaw = vang(F );

            if ~isnan(yaw) 
                obj.nodes(ch_id).yaw = yaw;
            end
        end
    end
    
    % adjust yaw according neighbors' weights
    for i=obj.sink+1:obj.n_nodes
        if obj.nodes(i).type ~= 'w'
            continue
        end
        force = [cos(obj.nodes(i).yaw), sin(obj.nodes(i).yaw)];
        %force = [0, 0];
        for j=obj.nodes(i).neighbor
            force = force + obj.nodes(j).weight*[cos(obj.nodes(j).yaw), sin(obj.nodes(j).yaw)];
        end
        yaw = atan(force(2)/force(1))+(force(1)<0)*pi;
        yaw = mod(yaw, 2*pi);
        if isnan(yaw)
            yaw = obj.nodes(i).yaw;
        end
        obj.nodes(i).yaw = yaw;
    end
    % Update cluster's direction for keeping connection
    for ch_id = obj.chs
        cpid = obj.nodes(ch_id).cpid;
        if cpid==0
            if isempty(obj.nodes(ch_id).force)
                obj.nodes(ch_id).v = obj.velocity;
            else
                FL = -obj.nodes(ch_id).force;
                yaw = obj.nodes(ch_id).yaw;
                vec = [cos(yaw), sin(yaw)];
                F = 10*vec+FL;
                obj.nodes(ch_id).yaw = vang(F);
                if obj.sink>0 && rand<0.1 && obj.trace{end}(ch_id).cpid==1
                    theta = vang(obj.nodes(1).p-obj.nodes(ch_id).p);
                    obj.nodes(ch_id).yaw = theta-pi/2+pi*rand();
                end
            end
        else
            if obj.sink>0 && cpid==1
                continue
            end
            yaw = obj.nodes(cpid).yaw;
            vec = [cos(yaw), sin(yaw)];
            %vecT = [cos(yaw), -sin(yaw)];
            vecT = [sin(yaw), -cos(yaw)];
            dist = vec * (obj.nodes(cpid).p-obj.nodes(ch_id).p);
            d = vecT * (obj.nodes(cpid).p-obj.nodes(ch_id).p);
            if dist<0
                obj.nodes(ch_id).v = obj.velocity;
            elseif dist>db
                obj.nodes(ch_id).v = obj.nodes(cpid).v+0.1;
            else
                obj.nodes(ch_id).v = obj.nodes(cpid).v;
            end
            a = 0.8;
            obj.nodes(ch_id).pd = obj.nodes(ch_id).pd*a + (abs(d)+1)*(1-a);
            FL = sign(d)*vecT;
            FL = FL*sign(abs(d)-obj.nodes(ch_id).pd);
            obj.nodes(ch_id).force = FL;
            F = 10*vec+FL;
            obj.nodes(ch_id).yaw = vang(F);
        end
    end
    
    % CMs follow the CH
    h = sqrt(max(obj.range_transmission^2 - 4*obj.range_coverage^2,0));
    for i=2:obj.n_nodes
        if obj.nodes(i).v==0
            continue
        end
        if obj.nodes(i).type == 'g'
            ch_id = obj.nodes(i).chid;
            %points = [obj.nodes([obj.nodes(ch_id).cluster_member, ch_id]).p];
            filter = [obj.nodes(obj.nodes(i).neighbor).chid]==ch_id;
            cm = obj.nodes(i).neighbor(filter);
            cm = [ch_id, i, cm];
            points = [obj.nodes(cm).p];

            % for original direction
            vec = [cos(obj.nodes(ch_id).yaw), sin(obj.nodes(ch_id).yaw)];
            dist = vec * points;
            [dist_ord, I] = sort(dist);
            i_ch = find(I==1, 1);
            i_self = find(I==2, 1);
            if i_self<i_ch % in the left of CH
                delta = dist_ord((i_self+1))-dist_ord((i_self));
                delta = max(0, delta-h);

            else    % in the right of CH
                delta = dist_ord((i_self))-dist_ord((i_self-1));
                delta = max(0, delta-h);
                delta = -delta;
            end
            v_P = obj.nodes(ch_id).v + 10*delta/(obj.range_transmission-h);

            % for orthogonal direction
            vec = [cos(obj.nodes(ch_id).yaw-pi/2), sin(obj.nodes(ch_id).yaw-pi/2)];
            dist = vec * points;
            [dist_ord, I] = sort(dist);
            i_ch = find(I==1, 1);
            i_self = find(I==2, 1);
            dd = min(obj.range_coverage*2, obj.range_transmission);
            if i_self<i_ch % in the left of CH
                delta = dist_ord((i_self+1))-dist_ord((i_self));
                delta = delta - dd;
            else    % in the right of CH
                delta = dist_ord((i_self))-dist_ord((i_self-1));
                delta = dd - delta;
            end
            v_T = 10*delta/(obj.range_coverage*2);

            obj.nodes(i).v = norm([v_P, v_T]);
            obj.nodes(i).yaw = mod(obj.nodes(ch_id).yaw - atan(v_T/v_P),2*pi);
        else
            obj.nodes(i).v = obj.velocity;
        end
        if isnan(obj.nodes(i).yaw)
            error('NAN yaw of Node-%d in following CH\n',i);
        end
    end
end

%% MPCA
function obj = clustering_MPCA(obj)
    % calculate the weights
    c = [1 1 1];
    for i=1:obj.n_nodes
        avg_LET = 0;
        avg_d = 0;
        for j=obj.nodes(i).neighbor
            avg_LET = avg_LET + link_exp_time(obj, i, j);
            avg_d = avg_d + length(obj.nodes(j).neighbor);
        end
        d = length(obj.nodes(i).neighbor);
        d = max(1, d);
        avg_LET = avg_LET / d;
        avg_d = avg_d / d;
        obj.nodes(i).weight_mpca = dot(c, [d, avg_d, log(1+avg_LET)]);
    end
    % 
    obj = clustering_mpca(obj);
    
    obj = join_nearest_cluster(obj);
    
    min_cm = 2;
    ch = [];
    for i=1:length(obj.chs)
        ch_id=obj.chs(i);
        if length(obj.nodes(ch_id).cluster_member)<min_cm...
                || (max([obj.nodes(obj.nodes(ch_id).cluster_member).weight_mpca])>obj.nodes(ch_id).weight_mpca && rand<0.1)
            for j=obj.nodes(ch_id).cluster_member
                obj.nodes(j) = obj.nodes(j).bewhite();
            end
            obj.nodes(ch_id) = obj.nodes(ch_id).bewhite();
            ch(end+1) = i;
        end
    end
    if ~isempty(ch)
        myprintf(' Delete CH: '); myprintf('%d ', obj.chs(ch)); myprintf('\n');
        obj.chs(ch) = [];
        myprintf(' CHs: '); myprintf('%d ', obj.chs); myprintf('\n');
    end
end

function obj = clustering_mpca(obj)
    for i=1:obj.n_nodes
        if obj.nodes(i).type == 'w'
            nb = obj.nodes(i).neighbor;
            wnb = nb([obj.nodes(nb).type]=='w');
            if obj.nodes(i).weight_mpca > max([obj.nodes(wnb).weight_mpca])
                obj.chs(end+1) = i;
                myprintf(' Add CH: %d\n', i);
                obj.nodes(i).type = 'b';
                obj.nodes(i).chid = i;
            end
        end
    end
    myprintf(' CHs: '); myprintf('%d ', obj.chs); myprintf('\n');
end

function t = link_exp_time(obj, i, j)
    % LET
    d = obj.nodes(j).p-obj.nodes(i).p;
    v = obj.nodes(j).v*[cos(obj.nodes(j).yaw); sin(obj.nodes(j).yaw)];
    v = v - obj.nodes(i).v*[cos(obj.nodes(i).yaw); sin(obj.nodes(i).yaw)];
    r = obj.range_transmission;
    t = (sqrt(dot(d,v)^2+dot(v,v)*(r^2-dot(d,d))))/dot(v,v);
    if t<0
        error('negtive')
    end
end

%% HDCA
function obj = clustering_HDCA(obj)
    
end

%% RWPMM
function obj = update_state_RWPMM(obj)
    pt = 1;
    wt = 10;
    pause_time = pt/obj.sim_time_step; 
    walk_time = wt/obj.sim_time_step;  
    p1 = 1/(pause_time+1);
    p2 = 1/(walk_time+1);
    for i=obj.sink+1:obj.n_nodes
        if obj.nodes(i).v == 0 && rand<p1 % pause
            obj.nodes(i).v = obj.velocity*(wt+pt)/wt;
            obj.nodes(i).yaw = rand * 2*pi;
        elseif obj.nodes(i).v>0 && rand<p2 % walk
            obj.nodes(i).v = 0;
        end
    end
end

%% CVBMM
function obj = update_state_CVBMM(obj)
    % Coverage-based mobility model
    % node 1 is the sink
    for i=obj.sink+1:obj.n_nodes
        nb = obj.nodes(i).neighbor;
        nb = nb(nb>1);
        if isempty(nb)
            continue
        end
        F = obj.nodes(i).p - [obj.nodes(nb).p];
        F = F ./ sum(F .* F);
        F = F+[cos(obj.nodes(i).yaw); sin(obj.nodes(i).yaw)]/obj.range_transmission;
        % update direction
        obj.nodes(i).yaw = vang(F);
    end
end

%% ABMM
function obj = update_state_ABMM(obj)
% alpha-based mobility model
    % update followship weights
    for i=1:obj.n_nodes
        nc = min(15, length(obj.nodes(i).neighbor));
        hc = min(99, obj.nodes(i).hop2sink);
        obj.nodes(i).weight = evalfis([nc, hc], obj.fis);
    end
    % adjust yaw according neighbors' weights
    yaw = [];
    for i=1:obj.n_nodes
        force = [cos(obj.nodes(i).yaw), sin(obj.nodes(i).yaw)];
        %force = [0, 0];
        for j=obj.nodes(i).neighbor
            force = force + obj.nodes(j).weight*[cos(obj.nodes(j).yaw), sin(obj.nodes(j).yaw)];
        end
        yaw(i) = atan(force(2)/force(1))+(force(1)<0)*pi;
        yaw(i) = mod(yaw(i), 2*pi);
        if isnan(yaw(i))
            yaw(i) = obj.nodes(i).yaw;
        end
        obj.nodes(i).yaw = yaw(i);
    end
end

%% CNBMM
function obj = update_state_CNBMM(obj)
    % node 1 is the sink
    for i=obj.sink+1:obj.n_nodes
        next_pi = move(obj, i);
        if obj.nodes(i).hop2sink < inf
            % if i is connected to the sink via one or more hops
            jid = obj.nodes(i).pid2sink;
            next_pj = move(obj, jid);
            if norm(next_pi-next_pj) > obj.range_transmission
                if jid==1 % if next hop is sink
                    theta = vang(obj.nodes(jid).p-obj.nodes(i).p);
                    obj.nodes(i).yaw = theta-pi/2+pi*rand();
                else
                    [~,v] = duration(obj, i, jid);
                    di = [cos(obj.nodes(i).yaw); sin(obj.nodes(i).yaw)]+v/norm(v);
                    obj.nodes(i).yaw = vang(di);
                end
            end
        else
            m = 0;
            for j=obj.nodes(i).neighbor
                next_pj = move(obj, j);
                if norm(next_pi-next_pj) <= obj.range_transmission
                    m = 0;
                    break
                end
                if m==0 || duration(obj, i, m)<duration(obj, i, j)
                    m = j;
                end
            end
            if m > 0
                [~,v] = duration(obj, i, m);
                di = [cos(obj.nodes(i).yaw); sin(obj.nodes(i).yaw)]+v/norm(v);
                obj.nodes(i).yaw = vang(di);
            end
        end
    end
 end
    
function [t, v] = duration(obj, i, j)
    % time when j leaves i
    d = obj.nodes(j).p-obj.nodes(i).p;
    v = obj.nodes(j).v*[cos(obj.nodes(j).yaw); sin(obj.nodes(j).yaw)];
    r = obj.range_transmission;
    t = (sqrt(dot(d,v)^2+dot(v,v)*(r^2-dot(d,d))))/dot(v,v);
    v = d+v*t;
    if t<0
        error('negtive')
    end
end

function theta = vang(vec)
    % calculate the angle of a vector 0~2Pi
    theta = atan(vec(2)/vec(1))+(vec(1)<0)*pi;
    theta = mod(theta, 2*pi);
end

function p = move(obj, i)
    t = obj.sim_time_step;
    p = obj.nodes(i).p + t * obj.nodes(i).v * [cos(obj.nodes(i).yaw); sin(obj.nodes(i).yaw)];
end

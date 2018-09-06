function [stat, obj] = statistics(obj)
    link2sink = 0;
    time2sink = 0;
    linkdur = 0;
    linknum = 0;
    chnum = 0;
    cmnum = 0;
    cmdur = 0;
    cdur = 0;
    whitenum = 0;
    raff = 0;
    conn = 0;

    nlink = 0;
    nch = 0;
    ncm = 0;
    nreaf = 0;

    adj_pre = zeros(obj.n_nodes);
    timesteps = length(obj.trace);
    step = 100;
    k=0;
    for t=1:timesteps
        nodes = obj.trace{t};
        adj_now = zeros(obj.n_nodes);

        % connectivity metrics
        % obj.connectivity(end+1) = sum([obj.nodes(1:obj.n_nodes).hop2sink]<inf);
        link2sink = link2sink + sum([nodes(1:obj.n_nodes).hop2sink]<inf);
        if mod(t,step)==0
            k=k+1;
            fprintf('%d ', k);
            %if  (strcmp(obj.model.mm, 'CNBMM')||strcmp(obj.model.mm, 'ABMM'))&&obj.n_nodes>100
                % do nothing
                %conn = conn + 1;
            %else
            con_t = connectivity(nodes);
            conn = conn + con_t;
            obj.connectivity(end+1) = con_t;
            %end
        end
        for i=1:obj.n_nodes
            linknum = linknum + length(nodes(i).neighbor);
            for j=nodes(i).neighbor
                if ~adj_pre(i,j)
                    nlink = nlink+1;
                end
                adj_now(i,j) = true;
            end

            % clustering metrics
            if nodes(i).type=='b'
                chnum = chnum + 1;
                if t==1 || obj.trace{t-1}(i).type~='b'
                    nch = nch+1;
                end
            elseif nodes(i).type=='g'
                cmnum = cmnum + 1;
                if t==1 || obj.trace{t-1}(i).type~='g'
                    ncm = ncm+1;
                end
            else
                whitenum = whitenum + 1;
            end
            %nodes(i).chid ~= obj.trace{t-1}(i).chid ||
            if t>1 && ( nodes(i).type ~= obj.trace{t-1}(i).type )
                nreaf = nreaf+1;
            end
        end
        adj_pre = adj_now;
    end

    time2sink = link2sink*obj.sim_time_step;
    linknum = linknum/2;
    linkdur = linknum*obj.sim_time_step;
    nlink = nlink/2;
    chdur = chnum*obj.sim_time_step;
    cmdur = cmnum*obj.sim_time_step;

    stat.cov_rate = obj.area_cover_rate(end);   % covarage rate
    stat.avgtime2sink = time2sink/obj.n_nodes;  % time link to sink
    stat.avglink2sink = link2sink/timesteps;    % # of nodes linked to sink
    stat.avglinkdur = linkdur/nlink;            % link duration
    stat.avgdegree = linknum/timesteps/obj.n_nodes; % average degree
    stat.avgconn = step*conn/timesteps;

    stat.avgchnum = chnum/timesteps;            % # of CHs
    stat.avgchdur = chdur/nch;                  % CHs duration
    stat.avgcmdur = cmdur/ncm;                  % CMs duration
    stat.avgorphannum = whitenum/timesteps;     % # of orphan/white nodes
    stat.avgreaff = nreaf/timesteps;            % Reaffiliation Frequency
    %stat.avgcmnum = cmnum/timesteps;            % # of CMs
    %stat.avgcmdur = cmdur/ncm;
end

function conn = connectivity(nodes)
    n = length(nodes);
    flag = zeros(1, n);
    cnt = [];
    conn = 0;
    k = 1;
    while k<=n
        if flag(k)
            k = k+1;
            continue
        end
        cnt(end+1) = 0;
        queue = {[k]};
        flag(k) = 1;
        k = k+1;
        while ~isempty(queue)
            cnt(end) = cnt(end)+length(queue{1});
            for i=queue{1}
                nb = nodes(i).neighbor;
                queue{end+1} = nb(~flag(nb));
                flag(nb) = 1;
            end
            queue(1)=[];
        end
    end
    cnt = cnt(cnt>1);
    if cnt
        conn = sum(arrayfun(@(x)nchoosek(x,2), cnt)) / nchoosek(length(nodes),2);
    end
end




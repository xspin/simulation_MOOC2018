function sys = fis_ABMM()
    sys = newfis('alpha');
    sys = addvar(sys, 'input', 'neighbor', [0, 15]);
    sys = addmf(sys,'input',1,'s(1)','gaussmf',[1 0]);
    sys = addmf(sys,'input',1,'m(2~3)','gaussmf',[1 2.5]);
    sys = addmf(sys,'input',1,'b(4~)','gaussmf',[5 15]);
    sys = addvar(sys, 'input', 'hop', [0, 99]);
    sys = addmf(sys,'input',2,'s(0)','gaussmf',[1 0]);
    sys = addmf(sys,'input',2,'m(1~2)','gaussmf',[1 2]);
    sys = addmf(sys,'input',2,'b(3~)','gaussmf',[3 10]);
    sys = addmf(sys,'input',2,'g(99)','gaussmf',[30 99]);
    sys = addvar(sys, 'output', 'alpha', [1, 10]);
    sys = addmf(sys, 'output',1,'s','gaussmf',[1 1]);
    sys = addmf(sys, 'output',1,'m','gaussmf',[1 5.5]);
    sys = addmf(sys, 'output',1,'b','gaussmf',[1 10]);
    % figure
    % plotmf(sys, 'input', 1)
    % figure
    % plotmf(sys, 'input', 2)
    % figure
    % plotmf(sys, 'output', 1)
    %    s m b g (hops)
    % s: b b m s
    % m: m b m s
    % b: s m s s
    % s m b g
    % 1 2 3 4
    [in1, in2] = meshgrid(1:3, 1:4);
    in1 = reshape(in1, [numel(in1), 1]);
    in2 = reshape(in2, [numel(in2), 1]);
    out = [3 3 2 1 2 3 2 1 1 2 1 1]';
    rules = [in1, in2, out, ones(length(out), 2)];
    sys = addrule(sys, rules);
    % showrule(sys)

    % plotfis(sys)
    % evalfis([1 1], sys)

end
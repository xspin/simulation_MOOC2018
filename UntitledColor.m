axesObjs = get(gcf, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
h = findobj(axesObjs, 'type', 'line');

flag = 'cover';
set(gcf, 'Position', [100 300 500 300]);
col = Color;
NameArray = {'Marker', 'Color', 'LineWidth', 'MarkerSize', 'MarkerFaceColor'};
ms = 7;
ls = 1;
if strcmp(flag, 'cluster')
    % for clustering
    % our approach, MP, HD
    ValueArray = {'d', col.red, ls, ms, col.red; ...
        '^', col.green, ls, ms, col.green; ...
        'o', col.blue, ls, ms, col.blue };
    set(h, NameArray, ValueArray);
    lg = {'Highest Degree', 'Mobility Prediction', 'MOOC'};
    %set(h(i), 'Marker', 'd', 'Color', col, 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', col)
elseif strcmp(flag, 'cover')
    % for cov&con
    % Our, AB, CNB, CVB, RWP
    ValueArray = {'d', col.red, ls, ms, col.red; ...
        '<', col.green, ls, ms, col.green; ...
        '>', col.blue, ls, ms, col.blue; ...
        's', col.purple, ls, ms, col.purple; ...
        'o', col.black, ls, ms, col.black};
    set(h, NameArray, ValueArray);
    lg = {'Random Waypoint', 'Coverage-based', ...
    'Connectivity-based', 'Alpha-based', 'MOOC'};
end
legend(lg)
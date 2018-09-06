classdef Color
    properties
        % 2f4ea1 292b78
        blue = [47, 78, 161]/255;
        dblue = [41, 43, 120]/255;
        % ee2123 7b1214
        red = [238, 33, 35]/255;
        dred = [123, 18, 20]/255;
        % 6cbd45 9cce6e
        green = [108, 189, 69]/255;
        dgreen = [156, 206, 110]/255;
        % b4509f
        purple = [180, 80, 159]/255;
        black = [110, 110, 110]/255;
% figure
% hold on
% x = 0:0.1:10;
% col = blue/255;
% plot(x, sin(x), 's-', 'Color', col, 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', col, 'MarkerIndices', 1:10:length(x))
% col = red/255;
% plot(x, cos(x), '^-', 'Color', col, 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', col, 'MarkerIndices', 1:10:length(x))
% col = green/255;
% plot(x, cos(x+1), 'o-', 'Color', col, 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', col, 'MarkerIndices', 1:10:length(x))
% col = purple/255;
% plot(x, sin(x+1), 'p-', 'Color', col, 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', col, 'MarkerIndices', 1:10:length(x))
% col = black/255;
% plot(x, sin(2*x), 'd-', 'Color', col, 'LineWidth', 2, 'MarkerSize', 12, 'MarkerFaceColor', col, 'MarkerIndices', 1:10:length(x))
    end
    methods
        function obj = Color()
        end
    end
end
function fig = pr_plot(x, y, title_text)
% generates filled staircase plot of probabilities over time period of x
%
% INPUTS:
% x - timetable object
% y - array of Prs of equal length to height(x)
% title_text is a string object

%%
% check if title is provided
if nargin < 3 || isempty(title_text)
    title_text = ""; % Default no title
end

% Create the stairstep plot
plot(x.Time, y,'blue');
ylabel("Probability");
xlabel("Time");
title(title_text);

% Remove top and right spines
ax = gca;
ax.Box = 'off';
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% Adjust y-axis limits
ylim([min(y), max(y)]);

% Get the x-vals and y-vals to fill plot
x_fill = x.Time([1, 1:end]);
y_fill = y([1, 1:end]);

% Fill the area underneath the plot with a solid colour
hold on;
area(x_fill, y_fill, 'FaceColor', 'blue', 'FaceAlpha', 0.3);
fig = gca;
hold off;

end
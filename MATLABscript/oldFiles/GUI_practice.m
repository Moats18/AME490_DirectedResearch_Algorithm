num = 1;

% Create a figure
f = figure();
% Create a panel to hold the plot axis
axisPanel = uipanel(f, 'Position', [0.2 0.2 0.3 1], 'BackgroundColor', [1 1 1]);

% Create a new axis on the panel
leftAxis = axes(axisPanel, 'Position', [0.2 0.2 0.7 0.7]);

% Create a different panel to hold the button and a second axis
buttonPanel = uipanel(f, 'Position', [0.5 0 0.5 1], 'BackgroundColor', [0.7 0.7 0.7]);
% Create an axis on the right panel


leftPlotButton = uicontrol(buttonPanel, 'Style', 'pushbutton', ...
    'String', 'Plot Left', 'Units', 'normalized', 'Position', [0.1 0 0.3 0.1], ...
    'Callback', @(src, event) buttonPressCallback(leftAxis, num));


function buttonPressCallback(thisAxis, num)
    % Plot on the specified axis
    plot(thisAxis, rand(1, num), rand(1, num), 'r.');
end
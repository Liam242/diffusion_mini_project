function slider_graph_example
    % Create a figure window
    fig = figure('Name', 'Interactive Sine Wave', ...
                 'NumberTitle', 'off', ...
                 'Position', [500 400 600 400]);

    % Create an axes for plotting
    ax = axes('Parent', fig, ...
              'Position', [0.1 0.3 0.85 0.65]);
    grid on
    hold on

    % Initial data
    x = linspace(0, 2*pi, 1000);
    freq = 1;
    y = sin(freq * x);

    % Plot initial sine wave
    plt = plot(ax, x, y, 'LineWidth', 2);
    title(ax, sprintf('Sine Wave with Frequency = %.2f Hz', freq))
    xlabel('x')
    ylabel('sin(f * x)')

    % Create slider control
    slider = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Min', 0.1, 'Max', 10, 'Value', freq, ...
        'Units', 'normalized', ...
        'Position', [0.1 0.1 0.8 0.05], ...
        'Callback', @(src, event) updatePlot(src, plt, ax, x));

    % Create a label for the slider
    uicontrol('Parent', fig, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.1 0.17 0.8 0.05], ...
        'String', 'Adjust Frequency (Hz)', ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
end

function updatePlot(slider, plt, ax, x)
    % Get the current slider value
    freq = get(slider, 'Value');
    % Update the y-data of the plot
    y = sin(freq * x);
    set(plt, 'YData', y);
    % Update the title
    title(ax, sprintf('Sine Wave with Frequency = %.2f Hz', freq))
end

slider_graph_example
% Code to animate refrigeration cycle from ssc_refrigeration
%% Plot Description:
%
% This figure shows the evolution of the fluid states in the refrigeration
% cycle over time. The four points on the refrigeration cycle (compressor
% inlet, condenser inlet, expansion valve inlet, and evaporator inlet) are
% plotted on the pressure-enthalpy diagram. The dotted contour lines
% indicates the temperature and the gray curve represents the saturation
% dome.

% Copyright 2013-2018 The MathWorks, Inc.

% Generate simulation results if they don't exist
if(~exist('simlog_ssc_refrigeration', 'var'))
    sim('ssc_refrigeration')
end

% Reuse figure if it exists, else create new figure
if ~exist('h3_ssc_refrigeration', 'var') || ...
        ~isgraphics(h3_ssc_refrigeration, 'figure')
    h3_ssc_refrigeration = figure('Name', 'ssc_refrigeration');
end
figure(h3_ssc_refrigeration)
clf(h3_ssc_refrigeration)

PlotPressureEnthalpyDiagram(simlog_ssc_refrigeration)




% Plot refrigeration cycle points from simulation result
function PlotPressureEnthalpyDiagram(simlog)

set(gcf, 'Name', 'ssc_refrigeration: Pressure-Enthalpy Diagram')

% Load fluid properties and saturation data for pressure-enthalpy diagram
load r134aPHDiagram.mat p_TLU h_TLU T_TLU p_sat h_sat

% Get simulation results
t = simlog.Compressor.Mass_Flow_Sensor.Mass_Energy_Flow_Rate_Sensor_2P.M.series.time;
p1 = simlog.Compressor.Cycle_Sensor_1.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
p2 = simlog.Compressor.Cycle_Sensor_2.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
p3 = simlog.Expansion_Valve.Cycle_Sensor_3.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
p4 = simlog.Expansion_Valve.Cycle_Sensor_4.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
h1 = simlog.Compressor.Cycle_Sensor_1.Thermodynamic_Properties_Sensor_2P.H.series.values('kJ/kg');
h2 = simlog.Compressor.Cycle_Sensor_2.Thermodynamic_Properties_Sensor_2P.H.series.values('kJ/kg');
h3 = simlog.Expansion_Valve.Cycle_Sensor_3.Thermodynamic_Properties_Sensor_2P.H.series.values('kJ/kg');
h4 = simlog.Expansion_Valve.Cycle_Sensor_4.Thermodynamic_Properties_Sensor_2P.H.series.values('kJ/kg');

% Interpolate results to obtain equal time steps
N = 500;
t_mov = linspace(t(1), t(end), N);
p_cycle = interp1(t, [p1 p2 p3 p4 p1], t_mov);
h_cycle = interp1(t, [h1 h2 h3 h4 h1], t_mov);

% Create pressure-enthalpy diagram with temperature contours
[contourC, contourh] = contour(h_TLU, p_TLU, T_TLU, 200:20:400, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':');
clabel(contourC, contourh, 'LabelSpacing', 432)
hold on

% Plot cycle points at t = 0
plot(h_sat, p_sat, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1)
handle_cycle = plot(h_cycle(1,:), p_cycle(1,:), 'color', [0 0.5 0], 'LineWidth', 1);
handle1 = plot(h_cycle(1,1), p_cycle(1,1), 'o', 'MarkerFaceColor', [0.7, 0.9, 1], 'MarkerEdgeColor', [0.7, 0.9, 1]);
handle2 = plot(h_cycle(1,2), p_cycle(1,2), '^', 'MarkerFaceColor', [1, 0.7, 0.7], 'MarkerEdgeColor', [1, 0.7, 0.7]);
handle3 = plot(h_cycle(1,3), p_cycle(1,3), 's', 'MarkerFaceColor', [1, 0.4 ,0.4], 'MarkerEdgeColor', [1, 0.4 ,0.4]);
handle4 = plot(h_cycle(1,4), p_cycle(1,4), 'd', 'MarkerFaceColor', [0.4, 0.4, 1], 'MarkerEdgeColor', [0.4 ,0.4, 1]);
hold off
legend([handle1 handle2 handle3 handle4], 'Point 1', 'Point 2', 'Point 3', 'Point 4', 'Location', 'southeast')
set(gca, 'YScale', 'log')
xlabel('Specific Enthalpy (kJ/kg)')
ylabel('Pressure (MPa)')
titleHandle = title(sprintf('Pressure-Enthalpy Diagram (t = %.f s)', 0));

% Create Play/Pause button for animation.
status = 'paused';
idxPaused = 1;
hButton = uicontrol('Style', 'pushbutton', 'String', 'Pause', ...
    'Units', 'normalized', 'Position', [0.13 0.94, 0.1, 0.05], ...
    'Callback', @(hObject, eventData) playAnimation);

% Play animation
playAnimation


    function playAnimation
        % Nested function to loop through time and set the cycle point data
        try
            if strcmp(status, 'playing')
                status = 'paused';
                set(hButton, 'String', 'Play')
                return
            end
            
            status = 'playing';
            set(hButton, 'String', 'Pause')
            
            % Plot pressure and specific enthalpy at cycle points.
            for i = idxPaused : length(t_mov)
                if strcmp(status, 'paused')
                    % Save state of the animation.
                    idxPaused = i;
                    return
                end
                set(handle_cycle, 'XData', h_cycle(i,:), 'YData', p_cycle(i,:))
                set(handle1, 'XData', h_cycle(i,1), 'YData', p_cycle(i,1))
                set(handle2, 'XData', h_cycle(i,2), 'YData', p_cycle(i,2))
                set(handle3, 'XData', h_cycle(i,3), 'YData', p_cycle(i,3))
                set(handle4, 'XData', h_cycle(i,4), 'YData', p_cycle(i,4))
                set(titleHandle, 'String', sprintf('Pressure-Enthalpy Diagram (t = %.f s)', t_mov(i)))
                drawnow
            end
            
            status = 'paused';
            set(hButton, 'String', 'Play')
            idxPaused = 1;
            
        catch ME
            % End gracefully if user closed figure during the animation.
            if ~strcmp(ME.identifier, 'MATLAB:class:InvalidHandle')
                rethrow(ME)
            end
        end
    end

end
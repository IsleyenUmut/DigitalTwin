% Code to plot vapor quality from ssc_refrigeration
%% Plot Description:
%
% This figure plots the vapor quality at each of the four points on the
% refrigeration cycle. It shows that when the compressor is turned on, the
% evaporator absorbs enough heat from the refrigerator compartment to
% completely vaporize the refrigerant. The condenser then brings the vapor
% quality down to about 0.02. Flash evaporation occurs in the expansion
% valve such that the refrigerant enters the evaporator at a vapor quality
% of about 0.4.

% Copyright 2013-2015 The MathWorks, Inc.

% Generate simulation results if they don't exist
if ~exist('simlog_ssc_refrigeration', 'var')
    sim('ssc_refrigeration')
end

% Reuse figure if it exists, else create new figure
if ~exist('h2_ssc_refrigeration', 'var') || ...
        ~isgraphics(h2_ssc_refrigeration, 'figure')
    h2_ssc_refrigeration = figure('Name', 'ssc_refrigeration');
end
figure(h2_ssc_refrigeration)
clf(h2_ssc_refrigeration)

plotQuality(simlog_ssc_refrigeration)



% Create plot from simulation results
function plotQuality(simlog)

set(gcf, 'Name', 'ssc_refrigeration: Refrigeration Cycle Vapor Qualities')

% Get simulation results
t = simlog.Compressor.Mass_Flow_Sensor.Mass_Energy_Flow_Rate_Sensor_2P.M.series.time;
x(:,1) = simlog.Compressor.Cycle_Sensor_1.Vapor_Quality_Sensor_2P.X.series.values;
x(:,2) = simlog.Compressor.Cycle_Sensor_2.Vapor_Quality_Sensor_2P.X.series.values;
x(:,3) = simlog.Expansion_Valve.Cycle_Sensor_3.Vapor_Quality_Sensor_2P.X.series.values;
x(:,4) = simlog.Expansion_Valve.Cycle_Sensor_4.Vapor_Quality_Sensor_2P.X.series.values;

% Plot results
handles(1) = subplot(2, 2, 1);
plot(t, x(:,3), 'Color', [1, 0.4, 0.4], 'LineWidth', 1)
grid on
ylabel('Vapor Quality')
title('Point 3: Valve Inlet')

handles(2) = subplot(2, 2, 2);
plot(t, x(:,2), 'Color', [1, 0.7, 0.7], 'LineWidth', 1)
grid on
ylabel('Vapor Quality')
title('Point 2: Condenser Inlet')

handles(3) = subplot(2, 2, 3);
plot(t, x(:,4), 'Color', [0.4, 0.4, 1], 'LineWidth', 1)
grid on
xlabel('Time (s)')
ylabel('Vapor Quality')
title('Point 4: Evaporator Inlet')

handles(4) = subplot(2, 2, 4);
plot(t, x(:,1), 'Color', [0.7, 0.9 ,1], 'LineWidth', 1)
grid on
xlabel('Time (s)')
ylabel('Vapor Quality')
title('Point 1: Compressor Inlet')

set(handles, 'YLim', [0 1])
linkaxes(handles)

end
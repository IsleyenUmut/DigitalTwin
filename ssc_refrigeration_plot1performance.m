% Code to plot simulation results from ssc_refrigeration
%% Plot Description:
%
% This figure plots the refrigeration cycle performance over time including
% the pressures, temperatures, energy flows, and mass flows. It shows that
% this refrigeration cycle operates at a compressor pressure ratio of about
% 5.5. The coefficient of performance, which is the ratio of the heat
% extracted to the compressor power input, is approximately 4.

% Copyright 2013-2017 The MathWorks, Inc.

% Generate simulation results if they don't exist
if ~exist('simlog_ssc_refrigeration', 'var')
    sim('ssc_refrigeration')
end

% Reuse figure if it exists, else create new figure
if ~exist('h1_ssc_refrigeration', 'var') || ...
        ~isgraphics(h1_ssc_refrigeration, 'figure')
    h1_ssc_refrigeration = figure('Name', 'ssc_refrigeration');
end
figure(h1_ssc_refrigeration)
clf(h1_ssc_refrigeration)

plotPerformance(simlog_ssc_refrigeration)



% Create plot from simulation results
function plotPerformance(simlog)

set(gcf, 'Position', [10 100 1000 500], 'Name', 'ssc_refrigeration: Cycle Performance')

% Get simulation results
t = simlog.Compressor.Mass_Flow_Sensor.Mass_Energy_Flow_Rate_Sensor_2P.M.series.time;
mdot(:,1) = simlog.Compressor.Mass_Flow_Sensor.Mass_Energy_Flow_Rate_Sensor_2P.M.series.values('g/s');
mdot(:,4) = simlog.Expansion_Valve.Mass_Flow_Sensor.Mass_Energy_Flow_Rate_Sensor_2P.M.series.values('g/s');
W = simlog.Compressor.Controlled_Mass_Flow_Rate_Source_2P.power.series.values('W');
Q = simlog.Refrigerator_Compartment.Heat_Flow_Sensor.Heat_Flow_Rate_Sensor.H.series.values('W');
p(:,1) = simlog.Compressor.Cycle_Sensor_1.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
p(:,2) = simlog.Compressor.Cycle_Sensor_2.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
p(:,3) = simlog.Expansion_Valve.Cycle_Sensor_3.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
p(:,4) = simlog.Expansion_Valve.Cycle_Sensor_4.Pressure_Internal_Energy_Sensor_2P.P.series.values('MPa');
Tc = simlog.Refrigerator_Compartment.Temperature_Sensor.Temperature_Sensor.T.series.values('K');
T(:,4) = simlog.Expansion_Valve.Cycle_Sensor_4.Thermodynamic_Properties_Sensor_2P.T.series.values('K');

% Plot results
handles(1) = subplot(2, 3, 1);
plot(t, p(:,1), 'Color', [0.7, 0.9, 1], 'LineWidth', 1)
hold on
plot(t, p(:,2), 'Color', [1, 0.7, 0.7], 'LineWidth', 1)
hold off
ylabel('Pressure (MPa)')
title('Compressor Pressures')
legend('p_{in}', 'p_{out}', 'Location', 'best')

handles(2) = subplot(2, 3, 2);
plot(t, -Q, 'Color', [0.7, 0.9, 1], 'LineWidth', 1)
ylabel('Heat Flow Rate (W)')
title('Compartment Heat Extracted')

handles(3) = subplot(2, 3, 3);
plot(t, mdot(:,1), 'Color', [1, 0.7, 0.7], 'LineWidth', 1)
hold on
plot(t, mdot(:,4), 'Color', [0.4 ,0.4 ,1], 'LineWidth', 1)
hold off
ylabel('Mass Flow Rate (g/s)')
title('Mass Flow Rate')
legend('Compressor', 'Valve', 'Location', 'best')

handles(4) = subplot(2, 3, 4);
plot(t, p(:,2)./p(:,1), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1)
xlabel('Time (s)')
ylabel('Pressure ratio')
title('Compressor Pressure Ratio')

handles(5) = subplot(2, 3, 5);
plot(t, W, 'Color', [1, 0.7, 0.7], 'LineWidth', 1)
xlabel('Time (s)')
ylabel('Power (W)')
title('Compressor Power')

handles(6) = subplot(2, 3, 6);
plot(t, T(:,4), 'Color', [0.4, 0.4, 1], 'LineWidth', 1)
hold on
plot(t, Tc, 'Color', [0.7, 0.9, 1], 'LineWidth', 1)
hold off
xlabel('Time (s)')
ylabel('Temperature (K)')
title('Evaporator Temperature')
legend('Inlet', 'Compartment', 'Location', 'best')

linkaxes(handles, 'x')

end
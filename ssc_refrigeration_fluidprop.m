function ssc_refrigeration_fluidprop(h1, h2)
% Code to plot fluid property data from ssc_refrigeration
%% Plot Description:
%
% The following two figures plot the fluid properties of refrigerant R-134a
% as a function of pressure (p) and specific internal energy (u) and as a
% function of pressure (p) and normalized internal energy (unorm),
% respectively. The fluid is a
%
% * subcooled liquid when -1 <= unorm < 0;
% * two-phase mixture when 0 <= unorm <= 1;
% * superheated vapor when 1 < unorm <= 2.
%
% The fluid property data is provided as a rectangular grid in p and unorm.
% Therefore, the grid in terms of p and u is non-rectangular.
%
% The R-134a fluid property data can be found in |r134aPropertyTables.mat|.

% Copyright 2013-2015 The MathWorks, Inc.

load r134aPropertyTables.mat r134aPropertyTables

r134aPropertyTables = interpolateTwoPhaseMixture(r134aPropertyTables);

plotNormalizedSpace(h1, r134aPropertyTables)
plotOriginalSpace(h2, r134aPropertyTables)

end


% Subfunction to interpolate fluid properties in the two-phase mixture
% region based on the liquid and vapor saturation properties
function fluidTables = interpolateTwoPhaseMixture(fluidTables)

% Create normalized internal energy vector for the two-phase mixture
mMixture = max(numel(fluidTables.liquid.unorm), numel(fluidTables.vapor.unorm));
fluidTables.mixture.unorm = linspace(0, 1, mMixture)';

% Linear interplation between liquid and vapor saturation
% Saturated liquid is the last row of the liquid property tables
% Saturated vapor is the first row of the vapor property tables
fluidTables.mixture.v  = interp1([0; 1], [fluidTables.liquid.v(end,:);  fluidTables.vapor.v(1,:)],  fluidTables.mixture.unorm);
fluidTables.mixture.s  = interp1([0; 1], [fluidTables.liquid.s(end,:);  fluidTables.vapor.s(1,:)],  fluidTables.mixture.unorm);
fluidTables.mixture.T  = interp1([0; 1], [fluidTables.liquid.T(end,:);  fluidTables.vapor.T(1,:)],  fluidTables.mixture.unorm);
fluidTables.mixture.nu = interp1([0; 1], [fluidTables.liquid.nu(end,:); fluidTables.vapor.nu(1,:)], fluidTables.mixture.unorm);
fluidTables.mixture.k  = interp1([0; 1], [fluidTables.liquid.k(end,:);  fluidTables.vapor.k(1,:)],  fluidTables.mixture.unorm);
fluidTables.mixture.Pr = interp1([0; 1], [fluidTables.liquid.Pr(end,:); fluidTables.vapor.Pr(1,:)], fluidTables.mixture.unorm);
fluidTables.mixture.u  = interp1([0; 1], [fluidTables.liquid.u_sat;     fluidTables.vapor.u_sat],   fluidTables.mixture.unorm);

end


% Subfunction to plot fluid properties as a function of pressure (p) and
% specific internal energy (u)
function plotOriginalSpace(hFigure, fluidTables)

% Concatenate liquid and vapor fluid property tables
u = [fluidTables.liquid.u; fluidTables.mixture.u; fluidTables.vapor.u];
p = repmat(fluidTables.p(:)', size(u, 1), 1);
fluidProps = {
    'Specific Volume (m^3/kg)',       [fluidTables.liquid.v;  fluidTables.mixture.v;  fluidTables.vapor.v ];
    'Specific entropy (kJ/kg/K)',     [fluidTables.liquid.s;  fluidTables.mixture.s;  fluidTables.vapor.s ];
    'Temperature (K)',                [fluidTables.liquid.T;  fluidTables.mixture.T;  fluidTables.vapor.T ];
    'Kinematic Viscosity (mm^2/s)',   [fluidTables.liquid.nu; fluidTables.mixture.nu; fluidTables.vapor.nu];
    'Thermal Conductivity (W/(m*K))', [fluidTables.liquid.k;  fluidTables.mixture.k;  fluidTables.vapor.k ];
    'Prandtl Number',                 [fluidTables.liquid.Pr; fluidTables.mixture.Pr; fluidTables.vapor.Pr];
    };

% Saturated liquid and vapor fluid properties
satProps = {
    fluidTables.liquid.v(end,:),  fluidTables.vapor.v(1,:);
    fluidTables.liquid.s(end,:),  fluidTables.vapor.s(1,:);
    fluidTables.liquid.T(end,:),  fluidTables.vapor.T(1,:);
    fluidTables.liquid.nu(end,:), fluidTables.vapor.nu(1,:);
    fluidTables.liquid.k(end,:),  fluidTables.vapor.k(1,:);
    fluidTables.liquid.Pr(end,:), fluidTables.vapor.Pr(1,:);
};

% Prepare figure
figure(hFigure)
clf(hFigure)
set(hFigure, 'Name', 'ssc_refrigeration: Fluid Properties', 'Toolbar', 'figure')
hAxes = gca;

% Create menu for selecting different fluid properties
hPopup = uicontrol('Style', 'popupmenu', 'String', fluidProps(:,1), ...
    'Units', 'normalized', 'Position', [0.331 0.932, 0.431, 0.05], ...
    'Value', 3, 'Callback', @(hObject, eventData) plotProperties, ...
    'FontWeight', 'bold');

% Plot fluid properties
plotProperties
view(hAxes, -37.5, 30)


    function plotProperties
        % Nested function to generate the surface plot for the selected
        % fluid property
        idx = get(hPopup, 'Value');
        [az, el] = view(hAxes);
        surf(hAxes, u, p, fluidProps{idx, 2})
        hold on
        plot3(hAxes, fluidTables.liquid.u_sat, fluidTables.p, satProps{idx,1}, 'k-', 'LineWidth', 2)
        plot3(hAxes, fluidTables.vapor.u_sat, fluidTables.p, satProps{idx,2}, 'k-', 'LineWidth', 2)
        hold off
        set(gca, 'YScale', 'log')
        view(hAxes, az, el)
        xlabel('Specific Internal Energy (kJ/kg)')
        ylabel('Pressure (MPa)')
        zlabel(fluidProps{idx, 1})
    end

end


% Subfunction to plot fluid properties as a function of pressure (p) and
% normalized internal energy (unorm)
function plotNormalizedSpace(hFigure, fluidTables)

% Concatenate liquid and vapor fluid property tables
n = numel(fluidTables.p);
p = fluidTables.p(:);
unorm = [fluidTables.liquid.unorm(:); fluidTables.mixture.unorm(:); fluidTables.vapor.unorm(:)];
[p, unorm] = meshgrid(p, unorm);
fluidProps = {
    'Specific Volume (m^3/kg)',       [fluidTables.liquid.v;  fluidTables.mixture.v;  fluidTables.vapor.v ];
    'Specific entropy (kJ/kg/K)',     [fluidTables.liquid.s;  fluidTables.mixture.s;  fluidTables.vapor.s ];
    'Temperature (K)',                [fluidTables.liquid.T;  fluidTables.mixture.T;  fluidTables.vapor.T ];
    'Kinematic Viscosity (mm^2/s)',   [fluidTables.liquid.nu; fluidTables.mixture.nu; fluidTables.vapor.nu];
    'Thermal Conductivity (W/(m*K))', [fluidTables.liquid.k;  fluidTables.mixture.k;  fluidTables.vapor.k ];
    'Prandtl Number',                 [fluidTables.liquid.Pr; fluidTables.mixture.Pr; fluidTables.vapor.Pr];
    };

% Saturated liquid and vapor fluid properties
satProps = {
    fluidTables.liquid.v(end,:),  fluidTables.vapor.v(1,:);
    fluidTables.liquid.s(end,:),  fluidTables.vapor.s(1,:);
    fluidTables.liquid.T(end,:),  fluidTables.vapor.T(1,:);
    fluidTables.liquid.nu(end,:), fluidTables.vapor.nu(1,:);
    fluidTables.liquid.k(end,:),  fluidTables.vapor.k(1,:);
    fluidTables.liquid.Pr(end,:), fluidTables.vapor.Pr(1,:);
};

% Prepare figure
figure(hFigure)
clf(hFigure)
set(hFigure, 'Name', 'ssc_refrigeration: Fluid Properties', 'Toolbar', 'figure')
hAxes = gca;

% Create menu for selecting different fluid properties
hPopup = uicontrol('Style', 'popupmenu', 'String', fluidProps(:,1), ...
    'Units', 'normalized', 'Position', [0.331 0.932, 0.431, 0.05], ...
    'Value', 3, 'Callback', @(hObject, eventData) plotProperties, ...
    'FontWeight', 'bold');

% Plot fluid properties
plotProperties
view(hAxes, -37.5, 30)


    function plotProperties
        % Nested function to generate the scatter plot for the selected
        % fluid property
        idx = get(hPopup, 'Value');
        [az, el] = view(hAxes);
        surf(hAxes, unorm, p, fluidProps{idx, 2})
        hold on
        plot3(hAxes, zeros(n, 1), fluidTables.p, satProps{idx,1}, 'k-', 'LineWidth', 2)
        plot3(hAxes, ones(n, 1), fluidTables.p, satProps{idx,2}, 'k-', 'LineWidth', 2)
        hold off
        set(gca, 'YScale', 'log')
        view(hAxes, az, el)
        xlabel('Normalized Internal Energy')
        ylabel('Pressure (MPa)')
        zlabel(fluidProps{idx, 1})
    end

end
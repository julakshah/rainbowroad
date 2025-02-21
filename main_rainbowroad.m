%% Main script for QEA2 rainbow road project
% AY '25 SPR
% Developers: Julian Shah | Nathaniel Banse
%% LEGEND
% VARIABLES
% u: symbolic variable in equations
% x: symbolic equation for i_hat position
% y: symbolic equation for j_hat position
% u_scope: the scope of u that the rainbow road has
% MAX_WHEEL_SPEED: the maximum wheel speed of the neatos
% WHEEL_BASE: the distance between neato wheels
% tu_vals: u values to plot tangents and normals at
% vx: symbolic x velocity equation
% vy: symbolic y velocity equation
% v_mag: symbolic velocity magnitude equation
% tx: symbolic normalized x tangent component equation
% ty: symbolic normalized y tangent component equation
% tx_vals: vector containing solved x tangent components
% ty_vals: vector containing solved y tangent components
% dtx: symbolic derivative of tx
% dty: symbolic derivative of ty
% dt_mag: symbolic tangent magnitude equation
% nx: symbolic normal vector x component
% ny: symbolic normal vector y component
% nx_vals: vector containing solved x normal components
% ny_vals: vector containing solved y normal components


% PLOTS

%% DEFINE CONSTANTS
% define parametric functions for rainbow road track
syms u
x = 0.3960 * cos(2.65 * (u + 1.4));
y = 0.99 * sin(u + 1.4);
% define the range of parametric curves
u_scope = [0, 3.2];

MAX_WHEEL_SPEED = .3; % m/s
WHEEL_BASE = .245; % meters

%% PLOT PATH/TANGENT/NORMAL

% TANGENT/NORMAL ORIGINS
% create "several points along the curve"
tu_vals = linspace(u_scope(1), u_scope(2), 5);
tnx = double(subs(x, u, tu_vals));
tny = double(subs(y, u, tu_vals));

% SOLVE FOR TANGENTS
% differentiate x and y functions
vx = diff(x);
vy = diff(y);
% find magnitude of v
v_mag = sqrt(vx^2 + vy^2);
% get x and y tangent equations
tx = vx / v_mag;
ty = vy / v_mag;
% substitute tu_vals into tx and ty
tx_vals = double(subs(tx, u, tu_vals));
ty_vals = double(subs(ty, u, tu_vals));

% SOLVE FOR NORMALS
% differentiate tx and ty
dtx = diff(tx);
dty = diff(ty);
% find magnitude of T
dt_mag = sqrt(dtx^2 + dty^2);
% get x and y normal equations
nx = dtx / dt_mag;
ny = dty / dt_mag;
% substitute tu_vals into nx and ny
nx_vals = double(subs(nx, u, tu_vals));
ny_vals = double(subs(ny, u, tu_vals));

%GRAPH
u_plotting = linspace(u_scope(1), u_scope(2));
x_vals = double(subs(x, u, u_plotting));
y_vals = double(subs(y, u, u_plotting));
figure(1); clf;
hold on
axis square
axis([-1.5,1.5,-1.5,1.5]);
title("Centerline Plot")
xlabel("x position")
ylabel("y position")
% plot the center line
plot(x_vals, y_vals, LineWidth=2);
% add dots to ends of center line
plot(x_vals(1),y_vals(1),'ko','markerfacecolor','k','markersize',5);
plot(x_vals(end),y_vals(end),'ko','markerfacecolor','k','markersize',5);
% plot tangent lines except the end points
for i = 2:length(tu_vals)-1
    plot(tnx(i),tny(i),'ko','markerfacecolor','k','markersize',5);
    quiver(tnx(i), tny(i), tx_vals(i), ty_vals(i), ...
        AutoScale="off",LineWidth=1.5, Color="r")

end
hold off

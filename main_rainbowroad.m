clear;
%% Main script for QEA2 rainbow road project
% AY '25 SPR
% Developers: Julian Shah | Nathaniel Banse
%% LEGEND

%{
VARIABLES
u: symbolic variable in equations
x: symbolic equation for i_hat position
y: symbolic equation for j_hat position
u_scope: the scope of u that the rainbow road has
MAX_WHEEL_SPEED: the maximum wheel speed of the neatos
WHEEL_BASE: the distance between neato wheels
tu_vals: u values to plot tangents and normals at
dx_du: symbolic x velocity 
dy_du: symbolic y velocity equation equation
v_mag: symbolic velocity magnitude equation
tx: symbolic normalized x tangent component equation
ty: symbolic normalized y tangent component equation
tx_vals: vector containing solved x tangent components
ty_vals: vector containing solved y tangent components
dtx: symbolic derivative of tx
dty: symbolic derivative of ty
dt_mag: symbolic tangent magnitude equation
nx: symbolic normal vector x component
ny: symbolic normal vector y component
nx_vals: vector containing solved x normal components
ny_vals: vector containing solved y normal components
T: Scalar that defines the total time we want driving the path to take
t: symbolic equation that defines t with respect to u
du_dt: scalar derived from t used to convert dx_du to velocity
v: symbolic equation converting to velocity with respect to time
d2x_du2: symbolic 2nd derivative of x
d2y_du2: symbolic 2nd derivative of y
cross_mag: symbolic equation for magnitude of curvature cross product
kappa: symbolic equation for curvature
omega: symbolic rotational velocity
vr: symbolic right wheel velocity
vl: symbolic left wheel velocity
t_vals: vector containing time for plotting
v_vals: vector containing velocity for plotting
omega_vals: vector containing angular velocity for graphing
vr_vals: vector containing right wheel velocity for plotting
vl_vals: vector containing left wheel velocity for plotting
%}


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
IP_STRING = '192.168.16.70';
NEATO_COMAPRE_FLAG = false;

%% PLOT PATH/TANGENT/NORMAL

% TANGENT/NORMAL ORIGINS
% create "several points along the curve"
tu_vals = linspace(u_scope(1), u_scope(2), 5);
tnx = double(subs(x, u, tu_vals));
tny = double(subs(y, u, tu_vals));

% SOLVE FOR TANGENTS
% differentiate x and y functions
dx_du = diff(x);
dy_du = diff(y);
% find magnitude of v
v_mag = sqrt(dx_du^2 + dy_du^2);
% get x and y tangent equations
tx = dx_du / v_mag;
ty = dy_du / v_mag;
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

% GRAPH
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
    quiver(tnx(i), tny(i), nx_vals(i), ny_vals(i), ...
        AutoScale="off",LineWidth=1.5, Color="g")

end
plot(tnx(1),tny(1),'ko','markerfacecolor','r','markersize',10);
hold off

%% DEFINE FUNCTIONS

% Time parameters
T = 30; % Total time (s), adjust if wheel speeds exceed 0.3 m/s
t = u * T / 3.2;
du_dt = 3.2 / T;

% Linear speed
v = v_mag * du_dt;

% Angular velocity via curvature
d2x_du2 = diff(dx_du);
d2y_du2 = diff(dy_du);
cross_mag = dx_du * d2y_du2 - dy_du * d2x_du2;
kappa = cross_mag / (v_mag^3);
omega = v .* kappa;

% Wheel velocities
vr = v + (omega * WHEEL_BASE) / 2;
vl = v - (omega * WHEEL_BASE) / 2;

%% PLOT PLANNED LINEAR SPEED AND ANGULAR VELOCITY VS. TIME

% Evaluating functions for plotting
t_vals = double(subs(t, u, u_plotting));
v_vals = double(subs(v, u, u_plotting));
omega_vals = double(subs(omega, u, u_plotting));
vr_vals = double(subs(vr, u, u_plotting));
vl_vals = double(subs(vl, u, u_plotting));

% Plot linear and angular velocities
figure('Position', [100, 100, 600, 400]);
subplot(2, 1, 1);
plot(t_vals, v_vals, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Linear Speed (m/s)');
title('Planned Linear Speed');
grid on;
subplot(2, 1, 2);
plot(t_vals, omega_vals, 'r-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
title('Planned Angular Velocity');
grid on;
sgtitle('Neatokart Center of Mass Velocities');
caption = sprintf('Planned linear speed and angular velocity of the Neatokart over %d s, derived from the Rainbow Road curve.', T);
text(0, -0.1, caption, 'FontSize', 8, 'HorizontalAlignment', 'center', 'Units', 'normalized');

% Plot wheel velocities
figure('Position', [700, 100, 600, 400]);
plot(t_vals, vr_vals, 'b-', 'LineWidth', 2, 'DisplayName', 'Right Wheel');
hold on;
plot(t_vals, vl_vals, 'r-', 'LineWidth', 2, 'DisplayName', 'Left Wheel');
plot([0 T], [MAX_WHEEL_SPEED MAX_WHEEL_SPEED], 'k--', 'DisplayName', 'Speed Limit');
plot([0 T], [-MAX_WHEEL_SPEED -MAX_WHEEL_SPEED], 'k--', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Wheel Velocity (m/s)');
title('Planned Wheel Velocities');
legend('Location', 'best');
grid on;
caption = sprintf('Left and right wheel velocities over %d s, constrained by ±0.3 m/s, for navigating Rainbow Road.', T);
text(T/2, -0.45, caption, 'FontSize', 8, 'HorizontalAlignment', 'center');
hold off;

%% NEATO MOVE

% connect to neato
neatov3.connect(IP_STRING);
neato_data = NaN; %this might have to be changed
% start the timer
tic;
% loop until 60 seconds have passed
while toc<=60
    % set t_in to the amount of time that has elapsed
    t_in = toc;
     
    % use MATLAB's linear interpolation to approximate
    % vl(t) and vr(t) at t=t_in
    vl_out = interp1(t_vals,vl_vals,t_in);
    vr_out = interp1(t_vals,vr_vals,t_in);
    if isnan(vl_out)
        vl_out = 0.0;
    end
    if isnan(vr_out)
        vr_out = 0.0;
    end
    neatov3.setVelocities(vr_out, vl_out)
    neato_data = [neato_data; neatov3.receive()];
end
neatov3.setVelocities(0.0, 0.0);
neato_data = neatov3.receive();
neatov3.disconnect()
%% NEATO COMPARE
% runs if NEATO_COMPARE


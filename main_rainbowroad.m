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
IP_STRING = '192.168.16.114';
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
axis equal
axis([-1,4,-1.5,1.5]);
title("Centerline Plot")
xlabel("X Position (m)")
ylabel("Y Position (m)")
% plot the center line
plot(x_vals, y_vals, LineWidth=2, DisplayName="Robot Centerline");
% add dots to ends of center line
plot(x_vals(1),y_vals(1),'ko','markerfacecolor','r','markersize',5, DisplayName="Path Start Point");
plot(x_vals(end),y_vals(end),'ko','markerfacecolor','k','markersize',5, DisplayName="Path End Point");
% plot tangent lines except the end points
for i = 2:length(tu_vals)-1
    quiver(tnx(i), tny(i), tx_vals(i), ty_vals(i), ...
        AutoScale="off",LineWidth=1.5, Color="r", DisplayName="Tangent Vector")
    quiver(tnx(i), tny(i), nx_vals(i), ny_vals(i), ...
        AutoScale="off",LineWidth=1.5, Color="g", DisplayName="Normal Vector")
end
legend()
caption = 'Planned path of Neato with tangent and normal vectors derived from Rainbow Road parametric equation';
text(1.3, -2.25, caption, 'FontSize', 8, 'HorizontalAlignment', 'center');
hold off

%% DEFINE FUNCTIONS

% Time parameters
T = 16; % Total time (s), adjust if wheel speeds exceed 0.3 m/s
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

% Evaluating functions for plotting
t_vals = double(subs(t, u, u_plotting));
v_vals = double(subs(v, u, u_plotting));
omega_vals = double(subs(omega, u, u_plotting));
vr_vals = double(subs(vr, u, u_plotting));
vl_vals = double(subs(vl, u, u_plotting));
%% NEATO MOVE

% connect to neato
neatov3.connect(IP_STRING);
neato_data = neatov3.receive();
t_in_data = 0;
% start the timer
tic;
% loop until 60 seconds have passed
while toc<=T
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
    t_in_data = [t_in_data; t_in]
    neato_data = [neato_data; neatov3.receive()];
end
neatov3.setVelocities(0.0, 0.0);
neatov3.disconnect()

%remove starting 0
t_in_data = t_in_data(2:end);
%% PLOT PLANNED LINEAR SPEED AND ANGULAR VELOCITY VS. TIME
neato_data_size = size(neato_data);
enc_data = zeros(neato_data_size(1), 2);
for i = 1:length(neato_data)
    enc_data(i, :) = neato_data(i).encoders();
end
left_wheel_encoder_list = enc_data(:, 2); % position data for left wheel
right_wheel_encoder_list = enc_data(:, 1); % position data for right wheel

dr = diff(right_wheel_encoder_list); % change in position for right wheel
dl = diff(left_wheel_encoder_list); % change in position for left wheel
theta = cumsum((dr - dl) ./ WHEEL_BASE); % theta vector
theta = theta + .316;
displacemennt = (dr + dl) ./ 2; % displacement vector

x_exp = cumsum(displacemennt .* cos(theta)); % x position vector
y_exp = cumsum(displacemennt .* sin(theta)); % y position vector

x_exp = x_exp - 0.33373270416764262492145011235346;
y_exp = y_exp + 0.97559523268857557885287983301804;


% SOLVE FOR TANGENTS
% differentiate x and y functions
dx_exp = diff(x_exp);
dy_exp = diff(y_exp);
% get x and y tangent equations
tx_exp = dx_exp ./ displacemennt(2:end, 1);
ty_exp = dy_exp ./ displacemennt(2:end, 1);

% SOLVE FOR NORMALS
% differentiate tx and ty
dtx_exp = diff(tx_exp);
dty_exp = diff(ty_exp);
% find magnitude of T
dt_mag_exp = sqrt(dtx_exp.^2 + dty_exp.^2);
% get x and y normal equations
nx_exp = dtx_exp ./ dt_mag_exp;
ny_exp = dty_exp ./ dt_mag_exp;

% Plot linear and angular velocities
figure('Position', [100, 100, 600, 400]); clf;
subplot(2, 1, 1);
hold on
plot(t_vals, v_vals, 'b-', 'LineWidth', 2, DisplayName='Planned Path');
plot(t_in_data(2:end), displacemennt(2:end)./diff(t_in_data), 'b--', 'LineWidth', 2, DisplayName='Experimental Path');
hold off
xlabel('Time (s)'); ylabel('Linear Speed (m/s)');
title('Planned Linear Speed');
legend()
grid on;
subplot(2, 1, 2);
hold on;
plot(t_vals, omega_vals, 'r-', 'LineWidth', 2, DisplayName='Planned Path');
plot(t_in_data(2:end), diff(theta)./diff(t_in_data), 'r--', 'LineWidth', 2, DisplayName='Experimental Path');
hold off;
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
title('Planned Angular Velocity');
legend()
grid on;
sgtitle('Neatokart Center of Mass Velocities');
caption = sprintf('Linear speed and angular velocity of the Neatokart over %ds, derived from the Rainbow Road curve.', T);
text(.5, -0.16, caption, 'FontSize', 8, 'HorizontalAlignment', 'center', 'Units', 'normalized');

% Plot wheel velocities
figure('Position', [700, 100, 600, 400]); clf;
hold on;
plot(t_vals, vr_vals, 'b-', 'LineWidth', 2, 'DisplayName', 'Right Wheel');
plot(t_vals, vl_vals, 'r-', 'LineWidth', 2, 'DisplayName', 'Left Wheel');
plot([0 T], [MAX_WHEEL_SPEED MAX_WHEEL_SPEED], 'k--', 'DisplayName', 'Speed Limit');
plot([0 T], [-MAX_WHEEL_SPEED -MAX_WHEEL_SPEED], 'k--', 'HandleVisibility', 'off');
xlabel('Time (s)'); ylabel('Wheel Velocity (m/s)');
title('Planned Wheel Velocities');
legend('Location', 'best');
grid on;
caption = sprintf('Left and right wheel velocities derived from angular and linear velocity over %ds, limited to ±0.3 m/s', T);
text(T/2, -0.25, caption, 'FontSize', 8, 'HorizontalAlignment', 'center');
hold off;

figure(); clf;
hold on
axis equal
axis([-1,4,-1.5,1.5]);
title("Experimental Vs. Theoretical Path")
xlabel("X position (m)")
ylabel("Y Position (m)")
% plot the center line
plot(x_vals, y_vals, LineWidth=2, DisplayName="Theoretical Path");
% add dots to ends of center line
plot(x_vals(1),y_vals(1),'ko','markerfacecolor','r','markersize',5, DisplayName="Path Start Point");
plot(x_vals(end),y_vals(end),'ko','markerfacecolor','k','markersize',5, DisplayName="Path End Point");
% plot tangent lines except the end points
for i = 2:length(tu_vals)-1
    quiver(tnx(i), tny(i), tx_vals(i), ty_vals(i), ...
        AutoScale="off",LineWidth=1.5, Color="#005000", DisplayName="Tangent Vector")
    quiver(tnx(i), tny(i), nx_vals(i), ny_vals(i), ...
        AutoScale="off",LineWidth=1.5, Color="#00cc00", DisplayName="Normal Vector")
end

plot(x_exp, y_exp, 'r--', LineWidth=2, DisplayName="Experimental Path")

for i = 1:3
    quiver(x_exp(20 * i + 15), y_exp(20 * i + 15), tx_exp(20 * i + 15), ty_exp(20 * i + 15), ...
        AutoScale="off",LineWidth=1.5, Color="#009900", LineStyle="--", DisplayName="Experimental Tangent Vector")
    quiver(x_exp(20 * i + 15), y_exp(20 * i + 15), nx_exp(20 * i + 15), ny_exp(20 * i + 15), ...
        AutoScale="off",LineWidth=1.5, Color="#cc9900", LineStyle="--", DisplayName="Experimental Normal Vector")
end
legend()
caption = 'Planned path of Neato compared to Measured path, with tangent and normal vectors at different points';
text(1.3, -2.25, caption, 'FontSize', 8, 'HorizontalAlignment', 'center');
hold off
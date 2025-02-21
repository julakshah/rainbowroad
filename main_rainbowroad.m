%% Main script for QEA2 rainbow road project
% AY '25 SPR
% Developers: Julian Shah | Nathaniel Banse

% Define parametric functions for rainbow road track
syms u
x = 0.3960 * cos(2.65 * (u + 1.4));
y = 0.99 * sin(u + 1.4);
% Define the range of parametric curves
u_scope = [0, 3.2];

MAX_WHEEL_SPEED = .3; % m/s
WHEEL_BASE = .245; % meters


%% PARAMETERS:
% number of agents
N = 3;
% dimension
d = 2;
% final time
T = 20;


%% SET THE RADIUS
R = N;
% hysteresis = 3*R/4;
Rh = 3*R/4;



%% INITIAL CONDITIONS
% initial positions
x0 = initx(N, d, N);
% x0 = x00; 

% initial velocities
v0 = initv(N, d, N);
% v0 = v00;




%% SET THE INCIDENCY MATRIX
global A;

A = set_adjacency_matrix(x0, N, R);











%% SOLVE THE PROBLEM FOR COMPARISON
arg0 = [reshape(x0', [N*d, 1]); reshape(v0', [N*d, 1])];

odefun = @(t, argx) F(argx, N, d, R, Rh);

[t, sol] = ode45(odefun, [0, T], arg0);





OUTPUT % script
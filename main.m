% clear all

%% PARAMETERS:
% number of agents
N = 3;
% dimension
d = 2;
% final time
T = 40;
% mesh length
n = 800;
% create mesch
mesh = Mesh(T, n);



%% INITIAL CONDITIONS
% initial positions
x0 = initx(N, d, N);
% x0 = x00;

% initial velocities
v0 = initv(N, d, N);
% v0 = v00;



%% SET OBJECTIVE PARAMETERS
alpha1 = 0; % integral of V(t)
alpha2 = 0; % V(T)
alpha3 = 0; % control
alpha4 = 0; % E(t)
alpha5 = 1; % integral of X(t)

%% CREATE THE DYNAMICS
gamma = 1;
delta = 1;
M = 1;
R = N;
dynamics = Dynamics(N, d, gamma, delta, alpha1, alpha3, alpha5, M, R);


%% CREATE THE OBJECTIVE
objective = Objective(dynamics, N, d, alpha2, alpha4);



%% CREATE RUNGE KUTTA SOLVER
A = [0 0 0; 0.5 0 0; -1 2 0];
b = [1.0/6.0    2.0/3.0    1.0/6.0];
% c = [0  0.5  1];
s = 3;
Nu = N;

arg0 = [reshape(x0', [N*d, 1]); reshape(v0', [N*d, 1]); 0];

rk = RungeKutta(A, b, s, dynamics, objective, arg0, 2*N*d+1, Nu, T, n);


%% INITIAL CONTROL GUESS
solu0 = ones(N, n,  s);


%% NCG MINIMIZATION
eps = 1;% not used 
sigma = 0.001;
limitLS = 5;
limitA = 25;
[solx, solu] = NCG(rk, objective, mesh, solu0, eps, sigma, limitLS, limitA);

sol = solx';
t = mesh.t;




%% SOLVE THE BFK PROBLEM FOR COMPARISON
soluBFK = ones(N, n,  s);
[solxBFK, solyBFK] = rk.solve_forward_equation(soluBFK);
solBFK = solxBFK';




OUTPUT; % script
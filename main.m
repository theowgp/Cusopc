% clear all

%% PARAMETERS:
% number of agents
N = 3;
% dimension
d = 2;
% final time
T = 50;
% mesh length
n = 100;
% create mesch
mesh = Mesh(T, n);



%% INITIAL CONDITIONS
% initial positions
x0 = initx(N, d, 4);

% initial velocities
v0 = initv(N, d, 1);







%% CREATE THE DYNAMICS
gamma = 1;
delta = 1;
M = 1;
R = 3;
dynamics = Dynamics(N, d, gamma, delta, M, R);



%% CREATE THE OBJECTIVE
objective = Objective(N, d);



%% CREATE RUNGE KUTTA SOLVER
A = [0 0 0; 0.5 0 0; -1 2 0];
b = [1.0/6.0    2.0/3.0    1.0/6.0];
% c = [0  0.5  1];
s = 3;
Nu = N;

arg0 = [reshape(x0', [N*d, 1]); reshape(v0', [N*d, 1]); 0];

rk = RungeKutta(A, b, s, dynamics, objective, arg0, 2*N*d+1, Nu, T, n);


%% INITIAL CONTROL GUESS
solu0 = zeros(N, n,  s);


%% NCG MINIMIZATION
eps = 1;% not used 
sigma = 0.001;
limitLS = 15;
limitA = 25;
[solx, solu] = NCG(rk, objective, mesh, solu0, eps, sigma, limitLS, limitA);

sol = solx';
t = mesh.t;




%% SOLVE THE BFK PROBLEM FOR COMPARISON
soluBFK = zeros(N, n,  s);
[solxBFK, solyBFK] = rk.solve_forward_equation(soluBFK);
solBFK = solxBFK';


% %% GET ENDTIME VALUES
% [xT, vT, zT, uT] = dynamics.convert(solx(:, end), solu(:, end, 1));

%% NORM of the SYSTEM VELOCITY at the end-time
normv = norm(solx(N*d+1:2*N*d, end))

%% NORM of the BFK SYSTEM VELOCITY at the end-time
normvBFK = norm(solxBFK(N*d+1:2*N*d, end))


%% PLOT THE LYAPUNOV FUNCTION
figure
for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YV(k) =  B(v, v, N);
end
plot(t, YV);
title('V(t) = 1/2N^2  sumij||vi -vj ||^2');


%% PLOT THE LYAPUNOV FUNCTION BFK
hold all
for k = 1:length(t)
%     x = reshape(solBFK(k, 1 : N*d), [d, N])';
    v = reshape(solBFK(k, N*d+1 : 2*N*d), [d, N])';
    YVBFK(k) =  B(v, v, N);
end
plot(t, YVBFK);
%  title('V(t) = 1/2N^2  sumij||vi -vj ||^2');





%% PLOT TRAJECTORIES
figure
for i = 1:N
    plot(sol(:, 2*i-1), sol(:, 2*i));
    hold all
end
title('evolution');

%% PLOT TRAJECTORIES BFK
figure
for i = 1:N
    plot(solBFK(:, 2*i-1), solBFK(:, 2*i));
    hold all
end
title('evolution BFK');



%% PLOT THE CONTROLS
% d = 1
figure
for i = 1:N
    plot(t(1:end-1), solu(i, :, 1));
    hold all
end
title(' norm of controls');


%% PLOT X
figure
for k = 1:length(t)
    x = reshape(sol(k, 1 : N*d), [d, N])';
    YX(k) =  B(x, x, N);
end
plot(t, YX);
title('X(t)');



%% PLOT E
figure
for k = 1:length(t)
    x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YE(k) =  E(x, v, N, R, dynamics);
end
plot(t, YE);
title('E(t)');
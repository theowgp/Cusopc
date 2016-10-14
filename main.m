% clear all

%% PARAMETERS:
% number of agents
N = 3;
% dimension
d = 2;
% final time
T = 100;
% mesh length
n = 300;
% create mesch
mesh = Mesh(T, n);



%% INITIAL CONDITIONS
% initial positions
x0 = initx(N, d, 10);
% x0 = x00; 

% initial velocities
v0 = initv(N, d, 2);
% v0 = v00;




%% CREATE THE DYNAMICS
gamma = 1;
delta = 1;
M = 1;
R = 3;
dynamics = Dynamics(N, d, gamma, delta, M, R);




%% CREATE RUNGE KUTTA SOLVER
A = [0 0 0; 0.5 0 0; -1 2 0];
b = [1.0/6.0    2.0/3.0    1.0/6.0];
% c = [0  0.5  1];
s = 3;
Nu = N;

arg0 = [reshape(x0', [N*d, 1]); reshape(v0', [N*d, 1])];

rk = RungeKutta(A, b, s, dynamics, arg0, 2*N*d, Nu, T, n);



%% SET THE TIME ARRAY
t = mesh.t;


%% SOLVE THE PROBLEM FOR COMPARISON
[solx, soly, solalpha] = rk.solve_forward_equation(1);
sol = solx';



%% SOLVE THE BFK PROBLEM FOR COMPARISON
[solxBFK, solyBFK, solalphaBFK] = rk.solve_forward_equation(0);
solBFK = solxBFK';


%% GET ENDTIME VALUES
[xT, vT] = convert(solx(:, end), N, d);


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




%% PLOT alphas
% d = 1
figure
for i = 1:N
    plot(t(1:end-1), solalpha(i, :, 1));
    hold all
end
title('alphas');



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
%% PLOT E BFK
hold all
for k = 1:length(t)
    x = reshape(solBFK(k, 1 : N*d), [d, N])';
    v = reshape(solBFK(k, N*d+1 : 2*N*d), [d, N])';
    YEBFK(k) =  E(x, v, N, R, dynamics);
end
plot(t, YEBFK);
title('E(t)');
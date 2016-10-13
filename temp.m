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
    plot(sol(1, 2*i-1), sol(1, 2*i), 'o');
    hold all
    plot(sol(:, 2*i-1), sol(:, 2*i));
    hold all
end
title('evolution');

%% PLOT TRAJECTORIES BFK
figure
for i = 1:N
    plot(solBFK(1, 2*i-1), solBFK(1, 2*i), 'o');
    hold all
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
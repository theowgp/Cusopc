%% CALCULATE ENERGIES
energy_my = control_energy(solx, dynamics, mesh, N, d, 'my')
energy_BFK = control_energy(solxBFK, dynamics, mesh, N, d, 'BFK')


%% PLOT THE LYAPUNOV FUNCTION
figure
for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YV(k) =  B(v, v, N);
end
plot(t, YV);
%% PLOT THE LYAPUNOV FUNCTION BFK
hold all
for k = 1:length(t)
%     x = reshape(solBFK(k, 1 : N*d), [d, N])';
    v = reshape(solBFK(k, N*d+1 : 2*N*d), [d, N])';
    YVBFK(k) =  B(v, v, N);
end
plot(t, YVBFK);
title('V(t)');





%% PLOT TRAJECTORIES
figure
ax = gca;
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
    plot(sol(1, 2*i-1), sol(1, 2*i), 'o');
    hold all
    plot(solBFK(:, 2*i-1), solBFK(:, 2*i));
    hold all
end
title('evolution BFK');




%% PLOT X
% figure
for k = 1:length(t)
    x = reshape(sol(k, 1 : N*d), [d, N])';
    YX(k) =  B(x, x, N);
end
plot(t, YX);
hold all
title('X(t)');




%% PLOT E
figure
for k = 1:length(t)
    x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YE(k) =  E(x, v, N);
end
plot(t, YE);
%% PLOT E BFK
hold all
for k = 1:length(t)
    x = reshape(solBFK(k, 1 : N*d), [d, N])';
    v = reshape(solBFK(k, N*d+1 : 2*N*d), [d, N])';
    YEBFK(k) =  E(x, v, N);
end
plot(t, YEBFK);
title('E(t)');


%% ANIMATED TRAJECTORIES PLOT
figure

% % draw initial conditions
% for i = 1:N
%     plot(sol(1, 2*i-1), sol(1, 2*i), 'o');
%     hold all
% end
% ax2 = gca;

%create animated lines
for i = 1:N
%     h(i) = animatedline(ax2, 'Color', 'b');
    h(i) = animatedline;
end
axis([ax.XLim ax.YLim]);


% start drawing animated lines
for k = 1:length(t)
    for i = 1:N
        addpoints(h(i), sol(k, 2*i-1), sol(k, 2*i));
        drawnow
    end
end
title('animated evolution');
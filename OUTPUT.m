
%% GET ENDTIME VALUES
[xT, vT] = convert(sol(end, :), N, d);
vT

%% NORM of the SYSTEM VELOCITY at the end-time
normv = norm(sol(end,  N*d+1:2*N*d))




%% PLOT THE LYAPUNOV FUNCTION
figure
for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YV(k) =  B(v, v, N);
end
plot(t, YV);
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




% %% PLOT X
% figure
% for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
%     YX(k) =  B(x, x, N);
% end
% plot(t, YX);
% title('X(t)');
% 
% 
% 
% %% PLOT E
% figure
% for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
%     v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
%     YE(k) =  E(x, v, N);
% end
% plot(t, YE);


% %% ANIMATED TRAJECTORIES PLOT
% figure
% 
% % draw initial conditions
% for i = 1:N
%     plot(sol(1, 2*i-1), sol(1, 2*i), 'o');
%     hold all
% end
% ax2 = gca;
% 
% %create animated lines
% for i = 1:N
%     h(i) = animatedline(ax2, 'Color', 'b');
% end
% axis([ax.XLim ax.YLim]);
% 
% 
% % start drawing animated lines
% for k = 1:length(t)
%     for i = 1:N
%         addpoints(h(i), sol(k, 2*i-1), sol(k, 2*i));
%         drawnow
%     end
% end
% title('animated evolution');
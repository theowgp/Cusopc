% %% PLOT THE LYAPUNOV FUNCTION
% figure
% for k = 1:length(t)
% %     x = reshape(sol(k, 1 : N*d), [d, N])';
%     v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
%     YV(k) =  B(v, v, N);
% end
% plot(t, YV);
% title('V(t) = 1/2N^2  sumij||vi -vj ||^2');
% 
% 
% %% PLOT THE LYAPUNOV FUNCTION BFK
% hold all
% for k = 1:length(t)
% %     x = reshape(solBFK(k, 1 : N*d), [d, N])';
%     v = reshape(solBFK(k, N*d+1 : 2*N*d), [d, N])';
%     YVBFK(k) =  B(v, v, N);
% end
% plot(t, YVBFK);
% %  title('V(t) = 1/2N^2  sumij||vi -vj ||^2');




% %% PLOT THE CONTROLS
% % d = 1
% figure
% for i = 1:N
%     plot(t(1:end-1), solu(2*i-1, :, 1));
%     hold all
% end
% title('controls d=1');
% % d = 2
% figure
% for i = 1:N
%     plot(t(1:end-1), solu(2*i, :, 1));
%     hold all
% end
% title('controls d=2');


% %% PLOT THE NORM CONTROLS
% % d = 1
% figure
% for i = 1:N
%     plot(t(1:end-1), solu(2*i-1, :, 1).^2 + solu(2*i, :, 1).^2);
%     hold all
% end
% title(' norm of controls');

% 
% 
%% PLOT E
figure
for k = 1:length(t)
    x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YE(k) =  E(x, v, N, R, dynamics);
end
plot(t, YE);
title('E(t)');

% 
% %% PLOT X
% figure
% for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
%     YX(k) =  B(x, x, N);
% end
% plot(t, YX);
% title('X(t)');

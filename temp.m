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


%% PLOT THE NORM CONTROLS
% d = 1
figure
for i = 1:N
    plot(t(1:end-1), solu(2*i-1, :, 1).^2 + solu(2*i, :, 1).^2);
    hold all
end
title(' norm of controls');

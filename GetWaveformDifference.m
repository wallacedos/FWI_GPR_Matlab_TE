function joint_source = GetWaveformDifference( isrc )
%CROSSCORRELATIONTIMESHIFTS Summary of this function goes here
%   Detailed explanation goes here

gather=0;
load(['Gather00_',num2str(isrc),'.mat'])
gather00 = gather;
figure()
imagesc( tout*1e9, recz, gather00');
ylabel('depth(m)')
xlabel('t(ns)')
title('Origin Data')
saveas(gcf,['gather00_',num2str(isrc),'.png'])


load(['Gather01_',num2str(isrc),'.mat'])
gather01 = gather;
figure()
imagesc(tout*1e9, recz, gather01');
ylabel('depth(m)')
xlabel('t(ns)')
title('Waveform Difference')
title('RTM Forward Data')
saveas(gcf,['gather01_',num2str(isrc),'.png'])


joint_source = gather00 - gather01;
% jointSource = zeros(size(gather_forward1));
% diff_gather2 = zeros(length(gather_forward1(:,1)),1);
% nsta = length(gather_forward1(1,:));
% tau = zeros(nsta,1);
% 
% std = (max(max(abs(gather_forward2(1:end-1,:)-gather_forward2(2:end,:))))/dt)^2 * 10;
% 
% for ista = 1:nsta
%     [dum, ndum] = max(abs(xcorr(gather_forward1(:,ista),gather_forward2(:,ista))));
%     tau(ista) = ndum - length(gather_forward1(:,1)) + 1;
%     diff_gather2(1:end-1) = (gather_forward2(1:end-1,ista) - gather_forward2(2:end,ista))/dt;
% %     if max(gather_forward2(:,ista)) < (max(max(gather_forward2))) * 0
% %         jointSource(:,ista) = tau(ista) * diff_gather2 * 0;
% %     else
% %         jointSource(:,ista) = tau(ista) * diff_gather2 / ((sum(diff_gather2.^2)));
% %     end
%     jointSource(:,ista) = tau(ista) * diff_gather2 / ((sum(diff_gather2.^2)) + std);
% end

figure();
imagesc( tout*1e9, recz, joint_source')
ylabel('depth(m)')
xlabel('t(ns)')
title('Waveform Difference')
joint_source = flipud(joint_source);
joint_source = joint_source';
saveas(gcf,['jointSource_',num2str(isrc),'.png'])

end


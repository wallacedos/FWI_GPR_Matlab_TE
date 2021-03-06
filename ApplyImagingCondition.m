function [ result ] = ApplyImagingCondition( f1, f2, isrc )

load(f1);
wavefield1 = wavefield;
load(f2);
wavefield2 = flipdim(wavefield, 1);

result = squeeze(sum(wavefield1 .* wavefield2, 1));

figure()
imagesc(x,z,result');
axis image


% nx = length(wavefield1(1,:,1));
% ny = length(wavefield1(1,1,:));
% 
% image1 = zeros(nx, ny);
% image2 = zeros(nx, ny);
% 
% 
% Ey11 = zeros(nx,ny);
% Ey21 = zeros(nx,ny);
% iEy = zeros(nx,ny);
% 
% for i  = 1:nx
%     for j = 1:ny
% %         iEy(i,j) = xcorr(Ey1(:,i,j),Ey2(end:-1:1,i,j),0);
%         iEy(i,j) = sum(Ey1(:,i,j) .* Ey2(end:-1:1,i,j));
%         
%     end
% %     disp(i);
% end
% figure()
% Ey_fig = iEy/max(max(abs(iEy)));
% % imagesc(iEy')
% xaxis = (1:nx)*dx;
% zaxis = (1:ny)*dz;
% 
% imagesc(xaxis,zaxis,iEy');
% % imagesc(xaxis,zaxis,Ey_fig');
% xlabel('m')
% ylabel('m')
% % caxis([-0.3,0.3])
% axis image;
% colorbar;
% saveas(gcf,['result_of_corr_',num2str(isrc),'.png'])
% 
% 
% figure()
% % imagesc(xaxis,zaxis,iEy');
% imagesc(xaxis,zaxis,Ey_fig');
% xlabel('m')
% ylabel('m')
% caxis([-0.3,0.3])
% axis image;
% % colorbar;
% saveas(gcf,['result_of_corr2_',num2str(isrc),'.png'])
% 
% 
% 
% save(['result_corr_',num2str(isrc),'.mat'],'iEy');
% 
% % nt = length(Ey(:,1,1));
% % for i = 1:3751
% %     figure(1)
% %     subplot(1,2,1)
% %     Ey11 = reshape(Ey1(i,:,:),nx,ny);
% %     imagesc(Ey11)
% %     subplot(1,2,2)
% %     Ey21 = reshape(Ey2(nt+1-i,:,:),nx,ny);
% %     imagesc(Ey21)
% %     pause(0.01);                                                                                                                                                                                                                                                              
% % end




end


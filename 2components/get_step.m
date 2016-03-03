function [ xistep, zistep ] = get_step( nsrc, inorm, iter)
%GET_STEP Summary of this function goes here
%   Detailed explanation goes here

xnumerator = 0;
xdenominator = 0;
znumerator = 0;
zdenominator = 0;
xresidual = 0;
zresidual = 0;
gather = 0;

for ii = 1:nsrc
%     ii = i*4;
% ii = fix(nsrc/2);
load(['Gather00_step_',num2str(ii),'.mat'])
xgather00_step = xgather;
zgather00_step = zgather;
% figure()
% imagesc(gather_forward1_step);
% saveas(gcf,['Gather_forward1_step_',num2str(ii),'.png'])

load(['Gather01_',num2str(ii),'.mat'])
xgather01 = xgather;
zgather01 = zgather;

%     gather_without_src = flipud(gather_forward1_step - gather_forward2);
xgather_without_src_step = xgather00_step - xgather01;
zgather_without_src_step = zgather00_step - zgather01;
% figure()
% imagesc(gather_without_src_step);
% saveas(gcf,['gather_without_src_step_',num2str(ii),'.png'])

load(['Gather00_',num2str(ii),'.mat'])
xgather00 = xgather;
zgather00 = zgather;

xgather_without_src = xgather00 - xgather01;
zgather_without_src = zgather00 - zgather01;
% figure()
% imagesc(gather_without_src);
% sum(sum(abs(gather_without_src))) / sum(sum(abs(gather_without_src_step)))
% figure()
% plot(gather_without_src(:,1))
% hold on
% plot(gather_without_src_step(:,1),'r')
% pause()
% 
% figure()
% plot(gather_forward1(:,1),'g')
% hold on
% plot(gather_forward1_step(:,1),'b')
% plot(gather_forward2(:,1),'r')
% legend('1','2','3')
% pause()
% sum(sum(abs(gather_without_src.*gather_without_src_step))) / sum(sum(abs(gather_without_src_step.^2)))

xdenominator = xdenominator + sum(sum(abs(xgather_without_src_step.^2)));
xnumerator = xnumerator + sum(sum((xgather_without_src.*xgather_without_src_step)));

zdenominator = zdenominator + sum(sum(abs(zgather_without_src_step.^2)));
znumerator = znumerator + sum(sum((zgather_without_src.*zgather_without_src_step)));

xresidual = xresidual +  sum(sum(abs(xgather_without_src.^2)));
zresidual = zresidual +  sum(sum(abs(zgather_without_src.^2)));

end

% istep = 0.1*(10/(iiii+10)) * sum(sum(abs(gather_without_src)))/sum(sum(abs(gather_without_src_step))) * norm;
% istep = 0.1 * sum(sum(abs(gather_without_src.*gather_without_src_step)))/sum(sum(abs(gather_without_src_step.^2))) * norm;
xistep =  xnumerator/xdenominator * inorm;
zistep =  znumerator/zdenominator * inorm;

fid=fopen('step.txt','a+');
fprintf(fid,'%d:: step: (%f, %f), numerator: (%f, %f), denominator: (%f, %f), residual: (%f, %f) \r\n',iter, xistep, zistep, xnumerator, xdenominator, znumerator, zdenominator, xresidual, zresidual);
fclose(fid)

end


for i = 1:size(u)
utemp = zeros(51,51);
utemp(1:51,1:51) = u(i,2:52,2:52);
%utemp(1,1) = 1;
imagesc(x(2:52),y(2:52),utemp)
set(gca,'YDir','normal')
colorbar
pause(.1)
end
%%
for i = 1:size(streamwise)
subplot(1,2,1)
plot(streamwise(i,:))
subplot(1,2,2)
plot(crosswise(i,:))
pause(0.01)
end
%%
utemp = zeros(51,51);
utemp(1:51,1:51) = u(58,2:52,2:52);
%utemp(1,1) = 1;
imagesc(x(2:52),y(2:52),utemp)
set(gca,'YDir','normal')
colorbar
title('Advection-Diffusion Final State')
%%
plot(mass)
function [] = waveplot(u,uexact)
% Plot movie of wave development over time
for i=1:size(u,1)
    plot(u(i,:))
    hold on
%    plot(uexact(i,:))
    hold off
    pause(0.01)
end
end
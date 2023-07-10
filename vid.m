video = VideoWriter('finalprojecttracer'); %open video file
video.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(video)

for i = 1:4:2*(Nt-1)/3+1
    zetatemp = zeros(Nx-2,Ny-2);
    zetatemp(1:Nx-2,1:Ny-2) = scalar1(i,2:Nx-1,2:Ny-1);
    pcolor(x(2:Nx-1),y(2:Ny-1),zetatemp')
    shading flat
    colorbar    
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(video, frame);
end
close(video)
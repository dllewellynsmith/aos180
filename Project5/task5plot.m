%%
figure(1)
imagesc(x(2:52),y(2:52),vorticity(2:52,2:52)')
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Vorticity Field')
%%
figure(2)
loglog(epsilon)
xlabel('Iteration')
ylabel('\epsilon')
title('Convergence')
%%
figure(3)
imagesc(x(2:52),y(2:52),phi(2:52,2:52)')
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Streamfunction')
%%
figure(4)
imagesc(x(2:52),y(2:52),u(2:52,2:52)')
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Horizontal Velocity Field')
%%
figure(5)
imagesc(x(2:52),y(2:52),v(2:52,2:52)')
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Vertical Velocity Field')
%%
figure(6)
imagesc(x(2:52),y(2:52),vort(2:52,2:52)')
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Laplacian of the Streamfunction')
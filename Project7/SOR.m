function [phi,count] = SOR(Nx,Ny,dx,dy,vorticity,tol,alpha)
num2 = max(max(abs(vorticity)));
num1 = 2/dx^2 + 2/dy^2;
phi = zeros(Nx,Ny);
R = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        R(i,j) = 1/dx^2 * (phi(i-1,j) - 2*phi(i,j)+phi(i+1,j)) + 1/dy^2 * (phi(i,j-1) - 2*phi(i,j)+phi(i,j+1)) + vorticity(i,j);
    end
end

epsilon = max(max(abs(R))) / ( num1*sum(abs(phi),'all') + num2);

count = 0;
while epsilon>tol
    count = count+1;
    for i = 2:Nx-1
            if i == Nx-1
                phi(Nx,:) = phi(3,:);
            end
        for j = 2:Ny-1
            R(i,j) = 1/dx^2 * (phi(i-1,j) - 2*phi(i,j)+phi(i+1,j)) + 1/dy^2 * (phi(i,j-1) - 2*phi(i,j)+phi(i,j+1)) + vorticity(i,j);
            phi(i,j) = phi(i,j) + alpha/num1 * R(i,j);
        end
    end
    phi(1,:) = phi(Nx-2,:);
    epsilon = max(max(abs(R))) / ( num1*sum(abs(phi),'all') + num2);
end
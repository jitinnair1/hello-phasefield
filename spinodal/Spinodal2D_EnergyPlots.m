clear all;
dt=0.5;
Nx=128;
Ny=128;
A=1.0;
plot_step=50;
noise=0.02;

more off;

% Declarations
conc=zeros(Nx,Ny);
energy1=zeros(1, plot_step);
energy2=zeros(1, plot_step);

% Initial profile
for i=1:Nx
    for j=1:Ny
        conc(i,j)=0.5 + noise*( 0.5 - rand() );
    end
end

mesh(conc);
view(2)
pause(0)

% Periodic Boundary
halfNx=Nx/2;
halfNy=Ny/2;
delkx=2*pi/Nx;
delky=2*pi/Ny;

% Evolve the profile


for z=1:plot_step
    for p=1:50
        
        g=2*A*conc.*(1-conc).*(1-2*conc);
        
        % FFT
        c_hat=fft2(conc);
        g_hat=fft2(g);
        
        
        for i=1:Nx
            for j=1:Ny
                
                %Periodic Boundary Condition
                
                %we take (i-1) to include the k = 0 point
                if ((i-1) <= halfNx)
                    kx=(i-1)*delkx;
                end
                
                if ((i-1) > halfNx)
                    kx=(i-1-Nx)*delkx;
                end
                
                if ((j-1) <= halfNy)
                    ky=(j-1)*delky;
                end
                
                if ((j-1) > halfNy)
                    ky=(j-1-Ny)*delky;
                end
                
                k2=kx*kx+ky*ky;
                k4=k2*k2;
                
                
                c_hat(i,j)=(c_hat(i,j)-dt*k2*g_hat(i,j))/(1+2*k4*dt);
                
            end
        end
        
        conc=real(ifft2(c_hat));
    end
    energy1(z)=cal_energy1(conc, Nx, Ny);
    energy2(z)=cal_energy2(conc, Nx, Ny);
%     mesh(conc);
%     view(2)
%     pause(0)
end

%% Variation of Bulk Energy
figure();
plot (energy1, 'r*');
title('Variation of Bulk Energy with successive iterations')
xlabel('Iteration No.'), ylabel('Bulk Energy')

%% Variation of Gradient Energy
figure();
plot (energy2, 'b*');
title('Variation of Gradient Energy with successive iterations')
xlabel('Iteration No.'), ylabel('Gradient Energy')
%print('Grad_Energy','-dpng')

%% Variation of Total Energy
figure();
total_energy=energy1+energy2;
plot (total_energy, 'g*');
title('Variation of Total Energy with successive iterations')
xlabel('Iteration No.'), ylabel('Total Energy')
%print('Total_Energy','-dpng')

%% Energy 1 Calculation
function [energy1] = cal_energy1(conc, Nx, Ny)
A=1.0;
energy1=0;

% energy1

for i=1:Nx
    for j=1:Ny
        energy1=energy1 + A*conc(i,j)*conc(i,j)*(1-conc(i,j))*(1-conc(i,j));
    end
end

end


%% Energy 2 Calculation
function [energy2] = cal_energy2(conc, Nx, Ny)
kappa=1.0;
energy2=0;

%energy 2

c_prime=zeros(Nx, Ny);

for i=1:Nx
    for j=1:Ny
        % Periodic condition for x
        wx=i-1;
        ex=i+1;
        if (wx<1)
            wx=wx+Nx;
        end
        if (ex>Nx)
            ex=ex-Nx;
        end
        
        % Periodic condition for y
        wy=j-1;
        ey=j+1;
        if (wy<1)
            wy=wy+Ny;
        end
        if (ey>Ny)
            ey=ey-Ny;
        end
        c_prime(i,j)=(conc(ex, j)-conc(wx, j))^2 +...
            (conc(i, ey)-conc(i, wy))^2;
        energy2 = energy2 + kappa*c_prime(i,j)*c_prime(i,j);
    end
end

end


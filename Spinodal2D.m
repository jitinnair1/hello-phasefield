clear all;
dt=0.5;
Nx=128;
Ny=128;
A=1.0;

more off;

% Declarations
conc=zeros(Nx,Ny);

% Initial profile
for i=1:Nx
    for j=1:Ny
        conc(i,j)=0.5 + ( 0.5 - rand() );
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


for z=1:50
    for p=1:20
        
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
    
    mesh(conc);
    view(2)
    pause(0)
end


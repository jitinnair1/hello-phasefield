clear all;
dt=0.1;
dx=1.0;
dy=1.0;
Nx=128;
Ny=128;
A=1.0;

more off;

% Declarations
conc=zeros(Nx,Ny);

% Initial profile
mflag=1; 

% Running 2D code in 1D mode to check for errors
% with random distribution only in x direction
if(mflag==0)
for i=1:1
    for j=1:Ny
        conc(i,j)=0.5 + ( 0.5 - rand() );
    end
end

for i=2:Nx
    for j=1:Ny
        conc(i,j)=conc(1,j);
    end
end
end

% Random distribution in both x and y direction
if(mflag==1)
for i=1:Nx
    for j=1:Ny
        conc(i,j)=0.5 + ( 0.5 - rand() );
    end
end
end

iprint=0; 
write_vtk_grid(Nx,Ny,dx,dy,iprint,conc);

% Periodic Boundary
halfNx=Nx/2;
halfNy=Ny/2;
delkx=2*pi/Nx;
delky=2*pi/Ny;

% Evolve the profile


for iprint=1:4
    for p=1:1000
        
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
    write_vtk_grid(Nx,Ny,dx,dy,iprint,conc);
end


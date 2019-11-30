% Details of how anisotropy is incorporated in the formulation
% can be found in DOI: 10.1080/01418610110038420

clear all;
dt=0.1;
dx=1.0;
dy=1.0;
Nx=128;
Ny=128;
A=1.0;
noise=0.02;
conc0=0.30;

%For cubic anisotropy
kappa=1.0;
gammaA=-20.0;
gammaI=25.0;

% Declarations
conc=zeros(Nx,Ny);

% Initial profile

% Random distribution in both x and y direction
for i=1:Nx
    for j=1:Ny
        conc(i,j)=conc0 + noise*( 0.5 - rand() );
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
                
                %including cubic anisotropy in interfacial energy
                kx4=kx*kx*kx*kx;
                ky4=ky*ky*ky*ky; 
                h_hat=kappa*k2 + gammaI*k4 + gammaA*(kx4 + ky4);
                
                %evolution equation including effects of anisotropy
                c_hat(i,j)=(c_hat(i,j)-dt*k2*g_hat(i,j))/(1+2*k4*dt*h_hat);
                
            end
        end
        
        conc=real(ifft2(c_hat));
    end
    write_vtk_grid(Nx,Ny,dx,dy,iprint,conc);
end


clear all;
dt=0.5;
Nx=64;
Ny=64;
Nz=64;
A=1.0;

more off;

% Declarations
conc=zeros(Nx,Ny, Nz);

% Initial profile
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
        conc(i,j,k) = 0.5 + ( 0.5 - rand() );
        end
    end
end


% Periodic Boundary
halfNx=Nx/2;
halfNy=Ny/2;
halfNz=Nz/2;
delkx=2*pi/Nx;
delky=2*pi/Ny;
delkz=2*pi/Nz;

% Evolve the profile


for z=1:2
    for p=1:5
        
        g=2*A*conc.*(1-conc).*(1-2*conc);
        
        % FFT
        c_hat=fftn(conc);
        g_hat=fftn(g);
        
        
        for i=1:Nx
            for j=1:Ny
                
                for b=1:Nz
                
                %Periodic Boundary Condition
                if ((i-1) <= halfNx) %we take (i-1) to include the k = 0 point
                    kx=(i-1)*delkx;
                end
                
                if ((i-1) > halfNx)
                    kx=(i-1-Nx)*delkx;
                end
                
                if ((j-1) <= halfNy) %we take (i-1) to include the k = 0 point
                    ky=(j-1)*delky;
                end
                
                if ((j-1) > halfNy)
                    ky=(j-1-Ny)*delky;
                end
                
                if ((b-1) <= halfNz) %we take (i-1) to include the k = 0 point
                    kz=(b-1)*delkz;
                end
                
                if ((b-1) > halfNz)
                    kz=(b-1-Nz)*delkz;
                end
                
                k2=kx*kx+ky*ky+kz*kz;
                k4=k2*k2;
                
                
                c_hat(i,j)=(c_hat(i,j)-dt*k2*g_hat(i,j))/(1+2*k4*dt);
                
                end
                
            end
        end
        
        conc=real(ifftn(c_hat));
    end
    
   
end


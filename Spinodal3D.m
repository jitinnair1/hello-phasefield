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

save_res(Nx, Ny, Nz, 0, conc);

% Periodic Boundary
halfNx=Nx/2;
halfNy=Ny/2;
halfNz=Nz/2;
delkx=2*pi/Nx;
delky=2*pi/Ny;
delkz=2*pi/Nz;

% Evolve the profile


for z=1:40
    for p=1:25
        
        g=2*A*conc.*(1-conc).*(1-2*conc);
        
        % FFT
        c_hat=fftn(conc);
        g_hat=fftn(g);
        
        
        for i=1:Nx
            for j=1:Ny
                for l=1:Nz
                
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
                
                if ((l-1) <= halfNz) %we take (i-1) to include the k = 0 point
                    kz=(l-1)*delkz;
                end
                
                if ((l-1) > halfNz)
                    kz=(l-1-Nz)*delkz;
                end
                
                k2=kx*kx+ky*ky+kz*kz;
                k4=k2*k2;
                
                
                c_hat(i,j,l)=(c_hat(i,j,l)-dt*k2*g_hat(i,j, l))/(1+2*k4*dt);
                
                end
                
            end
        end
        
        conc=real(ifftn(c_hat));
    end
    save_res(Nx, Ny, Nz, z, conc);
   
end

%% Save Results

function []=save_res(nx, ny, nz, istep, data1)
format long;

% Open Output File
fname=sprintf('time_%d.vtk', istep);
out=fopen(fname, 'w');
npoints=nx*ny*nz;

% Format header of VTK file
fprintf(out, '# vtk DataFile Version 2.0\n');
fprintf(out, 'time_10.vtk\n');
fprintf(out, 'ASCII\n');
fprintf(out, 'DATASET STRUCTURED_GRID\n');

% Co-ordinates of grid points
fprintf(out, 'DIMENSIONS %5d %5d %5d\n', nx, ny, nz);
fprintf(out, 'POINTS%7d float\n', npoints);

dx=1.0;
dy=1.0;
dz=1.0;

for i=1:nx
    for j=1:ny
        for k=1:nz
        x=(i-1)*dx;
        y=(j-1)*dy;
        z=(k-1)*dz;
        fprintf(out, '%14.6e %14.6e %14.6e\n', x, y, z);
        end
    end
end

% Write grid points
fprintf(out, 'POINT_DATA %5d\n', npoints);
fprintf(out, 'SCALARS CON float 1\n');
fprintf(out, 'LOOKUP_TABLE default\n');

for i=1:nx
    for j=1:ny
        for k=1:nz
        fprintf(out, '%14.6e\n', data1(i, j, k));
        end
    end
end
fclose(out);
end

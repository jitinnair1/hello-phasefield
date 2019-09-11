clear all;

% System Parameters
Nx=128;
Ny=128;
dx=0.5;
dy=0.5;
dt=0.1;
kappa=1.0;
L=1.0;
nstep=100;
nprint=10;
ttime=0.0;
dtime=0.005;


%Initial configuration
[etas, ngrain, glist] = initial_struc(Nx, Ny);

%Evolution
for istep=1:nstep
    ttime=ttime+dtime;
    for igrain=1:ngrain
        if (glist(igrain)==1)
            
            %Assign eta value to each grain
            for i=1:Nx
                for j=1:Ny
                    eta(i, j)=etas(i, j, igrain);
                end
            end
            
            %Evolve using periodic boundary
            for i=1:Nx
                for j=1:Ny
                    
                    %Periodic Boundary
                    iw=i-1;
                    ie=i+1;
                    jw=j-1;
                    je=j+1;
                    
                    if(iw<1)
                        iw=i-1+Nx;
                    end
                    
                    if(ie>Nx)
                        ie=i+1-Nx;
                    end
                    
                    if(jw<1)
                        jw=j-1+Ny;
                    end
                    
                    if(je>Ny)
                        je=j+1-Ny;
                    end
                    
                    %Calculate Laplacian
                    lap_eta(i,j) = (eta(ie,j) + eta(iw,j) + eta(i,je) + eta(i,jw)...
                        - 4.0*eta(i,j))/(dx*dy);
                    
                    %Calculate derivative of free energy
                    dfdeta=free_energy_diff(i, j, ngrain, etas, eta, igrain);
                    
                    %Integrate
                    eta(i,j)=eta(i,j)-dt*L*(dfdeta - kappa*lap_eta(i, j));
                    
                    %For small deviations
                    if (eta(i,j) >= 0.9999)
                        eta(i,j)=0.9999;
                    end
                    
                    if (eta(i,j) < 0.00001)
                        eta(i,j)=0.00001;
                    end
                    
                end
            end
            
            %Calculate Grain Volume
            grain_sum=0.0;
            
            for i=1:Nx
                for j=1:Ny
                    etas(i, j, igrain)=eta(i,j);
                    grain_sum=grain_sum + eta(i, j);
                end
            end
            
            %Check Grain Volume
            grain_sum=grain_sum/(Nx*Ny);
            if (grain_sum < 0.001)
                glist(igrain)=0;
            end
            
        end
    end
    
    %Print results every "nprint" steps
    if (mod(istep, nprint)==0)
        fname_area=sprintf('area_frac.txt');
        out2=fopen(fname_area, 'a');
        eta2=zeros(Nx, Ny);
        %fprintf(out2,'%14.6e',ttime);
        for igrain=1:ngrain
            ncount=0;
            for i=1:Nx
                for j=1:Ny
                    eta2(i, j)=eta2(i, j) + etas(i, j, igrain)^2;
                    if (etas(i, j, igrain) >= 0.5)
                        ncount=ncount+1;
                    end
                end
            end
            ncount=ncount/(Nx*Ny);
            fprintf(out2,'%14.6e\n', ncount);
        end
        save_res(Nx, Ny, dx, dy, istep, eta2);
    end
end

%% Initial Microstructure

function [etas, ngrain, glist] = initial_struc(Nx, Ny)
ngrain=2;

x0=Nx/2;
y0=Ny/2;
radius=20;

etas=zeros(Nx, Ny, ngrain);

for i=1:Nx
    for j=1:Ny
        etas(i, j, 1)=1.0;
        etas(i, j, 2)=0.0;
        if ((i-x0)*(i-x0)+(j-y0)*(j-y0)<radius*radius)
            etas(i, j, 1)=0.0;
            etas(i, j, 2)=1.0;
        end
    end
end

glist=zeros(1, ngrain);

for igrain=1:ngrain
    glist(igrain)=1.0;
end

end

%% Calculate derivative of free energy

function [dfdeta] = free_energy_diff(i, j, ngrain, etas, eta, igrain)
format long;
A=1.0;
B=1.0;
sum=0.0;

for jgrain=1:ngrain
    if (jgrain ~= igrain)
        sum=sum+etas(i, j, jgrain)*etas(i, j, jgrain);
    end
end

dfdeta=A*(2.0*B*eta(i, j)*sum + eta(i,j)^3 - eta(i, j));

end

%% Save Results

function []=save_res(nx, ny, dx, dy, istep, data1)
format long;

% Open Output File
fname=sprintf('time_%d.vtk', istep);
out=fopen(fname, 'w');
nz=1;
npoints=nx*ny*nz;

% Format header of VTK file
fprintf(out, '# vtk DataFile Version 2.0\n');
fprintf(out, 'time_10.vtk\n');
fprintf(out, 'ASCII\n');
fprintf(out, 'DATASET STRUCTURED_GRID\n');

% Co-ordinates of grid points
fprintf(out, 'DIMENSIONS %5d %5d %5d\n', nx, ny, nz);
fprintf(out, 'POINTS%7d float\n', npoints);

for i=1:nx
    for j=1:ny
        x=(i-1)*dx;
        y=(j-1)*dy;
        z=0.0;
        fprintf(out, '%14.6e %14.6e %14.6e\n', x, y, z);
    end
end

% Write grid points
fprintf(out, 'POINT_DATA %5d\n', npoints);
fprintf(out, 'SCALARS CON float 1\n');
fprintf(out, 'LOOKUP_TABLE default\n');

for i=1:nx
    for j=1:ny
        fprintf(out, '%14.6e\n', data1(i, j));
    end
end
fclose(out);
end

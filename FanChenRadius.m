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
ttime=0.0;
dtime=0.005;
r=zeros(1, nstep);


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
            
            r(istep)=get_radius(etas, Nx, Ny);
            
        end
    end
end

plot(r, '*')

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

%% Calculate Radius

function [radius] = get_radius(etas, Nx, Ny)
eta=zeros(Nx, Ny);
%Assign eta value to each grain
for i=1:Nx
    for j=1:Ny
        eta(i, j)=etas(i, j, 1);
    end
end
x0=Nx/2;
y0=Ny/2;
for m=1:y0
    if (eta(x0, m) < 1)
        rad_index_left=m;
        break
    end
end

for n=Ny:-1:y0
    if (eta(x0, n) < 1)
        rad_index_right=n;
        break
    end
end
radius = 0.5*(rad_index_right - rad_index_left);
end
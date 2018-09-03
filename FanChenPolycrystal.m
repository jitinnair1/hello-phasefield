clear all;

% System Parameters
Nx=64;
Ny=64;
dx=0.5;
dy=0.5;
dt=0.3;
kappa=1.0;
L=1.0;
nstep=1;
nprint=1;
ttime=0.0;
dtime=0.005;

%Declarartions
eta=zeros(Nx, Ny);
lap_eta=zeros(Nx, Ny);

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
                    eta(i,j)=eta(i,j)-dt*L*(dfdeta - L*kappa*lap_eta(i, j));

                    %For small deviations
                    if (eta(i,j) >= 0.9999)
                        eta(i,j)=0.999999999;
                    end

                    if (eta(i,j) < 0.00001)
                        eta(i,j)=0.000000001;
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
        eta2=zeros(Nx, Ny);
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
        end
        save_res(Nx, Ny, dx, dy, istep, eta2);
    end
end
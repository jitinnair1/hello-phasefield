profile on
clear sintering.m

% System Size Parametrs
Nx=200;
Ny=200;
dx=0.5;
dy=0.5;

% Time Integration Parametrs
nstep=50;
nprint=50;
dt=1.0e-4;

% Model Specific Parameters
coefm=5.0;
coefk=2.0;
coefl=5.0;
dvol=0.040;
dvap=0.002;
dsur=16.0;
dgrb=1.6;

% Prepare Microstructure
npart=2; % number of particless
eta=zeros(Nx, Ny); % for single active grain
phi=zeros(Nx, Ny);

% Initailse Laplacians
lap_eta=zeros(Nx, Ny);
lapcon1=zeros(Nx, Ny);
lapcon2=zeros(Nx, Ny);

% Get initial configuration
sflag=1;
[npart, etas, con]=micro_sint_pre(Nx, Ny, npart, sflag);


% Start Time Evolution
for istep=1:nstep

    % Sweep 1 for concentration using lapcon1
    for i=1:Nx
        for j=1:Ny

            ip = i-1;
            in = i+1;
            jp = j-1;
            jn = j+1;

            % First Laplacian Computation
            if (ip<1)
                ip=ip+Nx;
            end

            if (jp<1)
                jp=jp+Ny;
            end

            if (in>Nx)
                in=in-Nx;
            end

            if (jn>Ny)
                jn=jn-Ny;
            end

            lapcon1(i, j)=(con(ip, j) + con(in, j) + con(i, jn) + con(i, jp) - 4*con(i, j))/(dx*dy);

            % Free Energy Derivative for conc
            iflag=1; % for concentration derivative
            [dfcon, dfeta] = free_energy_diff(i, j, con, eta, etas, npart, iflag);

            % Second Laplacian Initailisation
            lapcon2(i, j) = dfcon - 0.5*coefm*lapcon1(i, j);
        end
    end %Nx Sweep1

    % Sweep 2 for concentration using lapcon2
    for i=1:Nx
        for j=1:Ny

            % Second Laplacian Computation
            ip = i-1;
            in = i+1;
            jp = j-1;
            jn = j+1;

            % First Laplacian Computation
            if (ip<1)
                ip=ip+Nx;
            end

            if (jp<1)
                jp=jp+Ny;
            end

            if (in>Nx)
                in=in-Nx;
            end

            if (jn>Ny)
                jn=jn-Ny;
            end

            lapcon2(i, j)=(con(ip, j) + con(in, j) + con(i, jn) + con(i, jp) - 4*con(i, j))/(dx*dy);

            % Order parameter relation
            phi=con(i, j)^3 * (10.0 - 15.0*con(i, j) + 6.0*con(i, j)^2);

            % Find total GB area
            sum=0.0;
            for ipart=1:npart
                for jpart=1:npart
                    if(ipart ~= jpart)
                        sum=sum+etas(i, j, ipart)*etas(i, j, jpart);
                    end % summation of grain boundary
                end
            end

            % Calculate Mobility
            mobil=dvol*phi + dvap*(1.0 - phi) + dsur*con(i, j)*(1-con(i,j)) + dgrb*sum;

            % Evolve conccentration
            con(i, j) = con(i, j) + dt*mobil*lapcon2(i, j);

            % For small devia(:, :)tions in concentration
            if con(i, j) >= 0.9999
                con(i, j)=0.9999;
            end

            if con(i, j) < 0.00001
                con(i, j)=0.00001;
            end

        end % Ny
    end % NX Sweep2

    % Sweep 3 for evolution of etas
    for ipart=1:npart
        eta=etas(:, :, ipart);

        for i=1:Nx
            for j=1:Ny

                % Laplacian for eta for ipart (lap_eta)
                ip = i-1;
                in = i+1;
                jp = j-1;
                jn = j+1;

                if (ip<1)
                    ip=ip+Nx;
                end

                if (jp<1)
                    jp=jp+Ny;
                end

                if (in>Nx)
                    in=in-Nx;
                end

                if (jn>Ny)
                    jn=jn-Ny;
                end

                lap_eta(i, j)=(eta(ip, j) + eta(in, j) + eta(i, jn) + eta(i, jp) - 4*eta(i, j))/(dx*dy);


                % Free Energy Derivative for etas
                iflag=2; % for eta derivative
                [dfcon, dfeta] = free_energy_diff(i, j, con, eta, etas, npart, iflag);

                % Evolve etas
                eta(i, j) = eta(i, j) - dt*coefl*(dfeta- 0.5 * coefk * lap_eta(i, j));

                % For small deviations in eta
                if eta(i, j) >= 0.9999
                    eta(i, j)=0.9999;
                end

                if eta(i, j) < 0.00001
                    eta(i, j)=0.00001;
                end
                
            end
        end

        %reassign eta to etas
        for i=1:Nx
            for j=1:Ny
                etas(i, j, ipart)=eta(i, j);
            end
        end
        
    end % ipart

    % Print Results
    if (mod(istep, nprint) == 0 || (istep == 1))
        fprintf("Done Step: %5d\n", istep);
        % write vtk file
        phi2=zeros(Nx, Ny);
        for i=1:Nx
            for j=1:Ny

                for ipart=1:npart
                    phi2(i, j) = phi2(i, j) + etas(i, j, ipart)^2;
                end
            end
        end %phi2 count

        write_vtk_grid_values(Nx, Ny, dx, dy, istep, con, phi2);

    end

end  % time iteration
profile off

% Save profile results for this run
profsave(profile('info'), 'Metrics/v4/')
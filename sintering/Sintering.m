clear all

% System Size Parametrs
Nx=256;
Ny=256;
dx=0.5;
dy=0.5;

% Time Integration Parametrs
nstep=5000;
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
iflag=1;
[npart, etas, con] = micro_sint_pre(Nx, Ny, npart, iflag);
eta=zeros(Nx,Ny);

% Calculate Laplacian
lap=delsq(numgrid('S', Nx)); %Assuming Nx=Ny

% Start Time Evolution
for istep=1:nstep
    
    % Sweep across Grid
    for i=1:Nx
        for j=1:Ny
            
            % Order parameter relation
            phi=con(i, j)^3 * (10.0 - 15.0*con(i, j) + 6.0*con(i, j)^2);
            
            % Find total GB area
            sum=0.0;
            for ipart=1:npart
                for jpart=1:npart
                    if(ipart ~= jpart)
                        sum=sum+etas(i, j, ipart)*etas(i, j, jpart);
                    end
                end
            end
            
            mobil=dvol*phi + dvap*(1.0 - phi) + dsur*con(i, j)*(1-con(i,j)) + drgb*sum;
            
            % Time Integration
            con(i, j) = con(i, j) + dt*mobil*lap(i, j);
            
            % For small deviations
            if con(i, j) >= 0.9999
                con(i, j)=0.9999;
            end
            
            if con(i, j) < 0.00001
                con(i, j)=0.00001;
            end
            
            
        end
    end
    
end

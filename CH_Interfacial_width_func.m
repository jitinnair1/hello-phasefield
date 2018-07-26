%This function is used to plot variation of IF width calculated with number
%of timesteps i.e. "nstep". Refer CH_width_variation for details.

%For this code solution of CH is done by FFT and calculation of slope is
%done by simple linear interpolation about the point of symmetry

function [width_spec] = CH_Interfacial_width_func(nstep)

dt=0.01;
N=400;
A=1.0;

% Declarations
conc=zeros(N,1);

% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        conc(i)=1;
    end
    if ((i==N/4) || (i==3*N/4))
        conc(i)=0.5;
    end
end


% Periodic Boundary
halfN=N/2;
delk=2*pi/N;

% Evolve the profile

for p=1:nstep
    
    % Define g
    g=2*A*conc.*(1-conc).*(1-2*conc);
    
    % FFT
    c_hat=fft(conc);
    g_hat=fft(g);
    
    for i=1:N
        
        %Periodic Boundary Condition
        if ((i) <= halfN)
            k=(i)*delk;
        end
        
        if ((i) > halfN)
            k=(i-N)*delk;
        end
        
        k2=k*k;
        k4=k2*k2;
        
        c_hat(i)=(c_hat(i)-dt*k2*g_hat(i))/(1+2*k4*dt);
    end
    
    conc=real(ifft(c_hat));
    
end

width_spec=get_N(conc, N);
end

function width = get_N(conc_vector, N)
lim=N/4;
lima=lim-1;
limb=lim+1;
m=(conc_vector(limb)-conc_vector(lima))/(limb-lima);
c0=conc_vector(lim)-(m*lim);
x2=(1.0-c0)/m;x1=-c0/m;
width=x2-x1;

end

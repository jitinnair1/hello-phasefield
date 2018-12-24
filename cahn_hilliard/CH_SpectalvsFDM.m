%% Spectral
clear all;
D=1.0;
dt=0.1;
dx=1;
N=128;
m=2;
A=1.0;
kappa=1.0;
energy1=0;
energy2=0;
energy3=0;
energy4=0;
nstep=1000;
beta1=dt/dx*dx;
beta2=2*kappa*beta1/dx*dx;
col_labels={'bulk';'gradient';'total'}; 

% Declarations
concSpectral=zeros(N,1);
c_prime=zeros(N, 1);

% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        concSpectral(i)=1;
    end
end

% Define g
g=zeros(N,1);


% Periodic Boundary
halfN=N/2;
delk=2*pi/N;

% Evolve the profile

for p=1:nstep
    
    % Define g
    g=2*A*concSpectral.*(1-concSpectral).*(1-2*concSpectral);
    
    % FFT
    c_hat=fft(concSpectral);
    g_hat=fft(g);
    
    for i=1:N
        
        %Periodic Boundary Condition
        if ((i-1) <= halfN)
            k=(i-1)*delk;
        end
        
        if ((i-1) > halfN)
            k=(i-1-N)*delk;
        end
        
        k2=k*k;
        k4=k2*k2;
        
        c_hat(i)=(c_hat(i)-dt*k2*g_hat(i))/(1+2*k4*dt);
        
        
        concSpectral=real(ifft(c_hat));
    end    
end

% Calculation of Spectral Interfacial energy

% energy1

for i=1:N
    energy1=energy1 + A*concSpectral(i)*concSpectral(i)...
        *(1-concSpectral(i))*(1-concSpectral(i));
end

% energy 2

%Periodic Boundary Condition
c_hat=fft(concSpectral);
for i=1:N

    if ((i-1) <= halfN)
        k=(i-1)*delk;
    end
    
    if ((i-1) > halfN)
        k=(i-1-N)*delk;
    end
    
    % Transform to fourier space for calculating derivative
    c_hat(i)=c_hat(i)*complex(0,1)*k;  
end

c_prime=real(ifft(c_hat));

for i=1:N
    energy2 = energy2 + kappa*c_prime(i)*c_prime(i);
end

E1S=0.5*energy1;
E2S=0.5*energy2;
E3S=0.5*(energy1 + energy2);
EnergySpectral=[E1S; E2S; E3S];

%% FDM

% Declarations
concFDM=zeros(N,1);
conc_old=zeros(N,1);

% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        conc_old(i)=1;
    end
end

% Define g
g=zeros(N,1);

% Evolve the profile


for k=1:nstep
    
    for i=1:N
        
        % Define g
        g=2*A*conc_old.*(1-conc_old).*(1-2*conc_old);
        
        w=i-1;
        ww=i-2;
        e=i+1;
        ee=i+2;
        
        if (ww<1)
            ww=ww+N;
        end
        if (w<1)
            w=w+N;
        end
        if (ee>N)
            ee=ee-N;
        end
        if(e>N)
            e=e-N;
        end
        
        % ellipsis are for text wrap in MATLAB editor window
        
        concFDM(i)=conc_old(i) + beta1*(g(w)-2*g(i)+g(e))...
            - beta2*(conc_old(ww)-4*conc_old(w)+6*conc_old(i)...
            -4*conc_old(e) + conc_old(ee));
        
        conc_old(i)=concFDM(i);
        
    end
    
end
plot(conc_old, 'b*');
hold on
plot(concSpectral, 'r*');
hold on

% Calculation of Interfacial energy

for i=1:N
    energy3=energy3 + A*conc_old(i)*conc_old(i)*(1-conc_old(i))*(1-conc_old(i));
end

% Interfacial For FDM
for i=1:N
    w=i-1;
    e=i+1;
    if (w<1)
        w=w+N;
    end
    if (e>N)
        e=e-N;
    end
    c_prime(i)=(conc_old(e)-conc_old(w))/(2*dx);
    energy4 = energy4 + kappa*c_prime(i)*c_prime(i);
end


E1F=0.5*energy3;
E2F=0.5*energy4;
E3F=0.5*(energy3 + energy4);
EnergyFDM=[E1F; E2F; E3F];


%% Plot of Delta
difference=zeros(1, N);

for i=1:N
    difference(i)=concSpectral(i)-conc_old(i);
end
clf;
plot(difference, 'g*')

%% Interfacial Energy - FDM vs Spectral

table(col_labels, EnergyFDM, EnergySpectral)




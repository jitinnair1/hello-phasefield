clear all;
D=1.0;
dt=0.5;
N=128;
m=2;
A=1.0;
kappa=1.0;

% Declarations
conc=zeros(N,1);
c_prime=zeros(N, 1);

% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        conc(i)=1;
    end
end

plot(conc, 'r')
hold on

% Define g
g=zeros(N,1);


% Periodic Boundary
halfN=N/2;
delk=2*pi/N;

% Evolve the profile

for p=1:60
    
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
plot(conc);

% Calculation of Interfacial energy

energy1=0;
energy2=0;

% energy1

for i=1:N
    energy1=energy1 + A*conc(i)*conc(i)*(1-conc(i)*conc(i))*(1-conc(i)*conc(i));
end

% energy 2

c_hat=fft(conc);

%Periodic Boundary Condition
for i=1:N
    
    if ((i) <= halfN)
        k=(i)*delk;
    end
    
    if ((i) > halfN)
        k=(i)*delk;
    end
    
    % Transform to fourier space for calculating derivative
    c_prime(i)=real(ifft(c_hat(i)*complex(0,1)*k));
    
end

for i=1:N
    energy2 = energy2 + kappa*c_prime(i)*c_prime(i);
end

E1S=0.5*energy1;
E2S=0.5*energy2;
E3S=0.5*(energy1 + energy2);

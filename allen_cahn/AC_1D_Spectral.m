clear all;
D=1.0;
dt=0.1;
N=128;
m=2;
L=1.0;
A=1.0;
kappa=1.0;

% Declarations
phi=zeros(N,1);

% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        phi(i)=1;
    end
end

plot(phi, 'r')
hold on

% Define g
g=zeros(N,1);

% Periodic Boundary
halfN=N/2;
delk=2*pi/N;

% Evolve the profile

for p=1:5
    
    for i=1:N
        
        %Define g
        g=2*A*phi.*(1-phi).*(1-2*phi);
        
        % FFT
        phi_hat=fft(phi);
        g_hat=fft(g);
        
        %Periodic Boundary Condition
        if ((i) <= halfN) %we take (i-1) to include the k = 0 point
            k=(i)*delk;
        end
        
        if ((i) > halfN)
            k=(i-N)*delk;
        end
        
        k2=k*k;
        
        phi_hat(i)=(phi_hat(i)-L*dt*g_hat(i))/(1+2*kappa*L*k2*dt);
        
        phi=real(ifft(phi_hat));
    end
    
    
end
plot(phi);


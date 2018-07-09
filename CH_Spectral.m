clear all;
D=1.0;
dt=0.1;
N=128;
m=2;
A=1.0;

% Declarations
conc=zeros(N,1);

% Initial profile
for i=1:N
    conc(i)=0.5*(1+sin(2*pi*m*i/N));
    %conc(i)=sin(2*pi*m*i/N);
end

plot(conc, 'r*')
hold on

% Define g
g=zeros(N,1);

% Periodic Boundary
halfN=N/2;
delk=2*pi/N;

% Evolve the profile

for p=1:3000
    
    % Define g
    g=2*A*conc.*(1-conc).*(1-2*conc);
    
    % FFT
    c_hat=fft(conc);
    g_hat=fft(g);
    
    for i=1:N
        
        %Periodic Boundary Condition
        if ((i-1) <= halfN) %we take (i-1) to include the k = 0 point
            k=(i)*delk;
        end
        
        if ((i-1) > halfN)
            k=(i-1-N)*delk;
        end
        
        k2=k*k;
        k4=k2*k2;
        
        c_hat(i)=(c_hat(i)-dt*k2*g_hat(i))/(1+2*k4*dt);
    end
    
    conc=real(ifft(c_hat));
end
plot(conc);


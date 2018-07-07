clear all;
D=1.0;
dt=0.001;
N=128;
m=1;
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

for i=1:N
    g(i)=2*A*conc(i).*(1-conc(i)).*(1-2*conc(i));
end

% Periodic Boundary
halfN=N/2;
delk=2*pi/N;

% Evolve the profile

% FFT
c_hat=fft(conc);
g_hat=fft(g);

for j=1:2
    for p=1:10000
        
        for i=1:N
            
            %Periodic Boundary Condition
            if ((i) <= halfN) %we take (i-1) to include the k = 0 point
                k=(i)*delk;
            end
            
            if ((i) > halfN)
                k=(i)*delk;
            end
            
            k2=k*k;
            k4=k2*k2;
            
            c_hat(i)=(c_hat(i)-dt*k2*g_hat(i))/(1+2*k4*dt);
        end
        conc=real(ifft(c_hat));
    end
    plot(conc);
end

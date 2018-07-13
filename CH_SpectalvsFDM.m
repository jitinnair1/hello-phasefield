clear all;
D=1.0;
dt=0.5;
N=64;
m=2;
A=1.0;
kappa=1.0;

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

for p=1:60
    
    % Define g
    g=2*A*concSpectral.*(1-concSpectral).*(1-2*concSpectral);
    
    % FFT
    c_hat=fft(concSpectral);
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
        
        
        concSpectral=real(ifft(c_hat));
    end
    
    
    
end
subplot(2,1,1)
plot(concSpectral, 'r*');

hold on

%% FDM

D=1.0;
dx=0.1;
dt=0.01;
m=2;
kappa=1.0;
A=1.0;
beta1=dt/dx*dx;
beta2=2*kappa*beta1/dx*dx;

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


for k=1:8000
    
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
subplot(2,1,1)
plot(conc_old, 'b*');


%% Plot of Delta
difference=zeros(1, N);

for i=1:N
    difference(i)=concSpectral(i)-conc_old(i);
end
subplot(2,1,2)
plot(difference, 'g*')



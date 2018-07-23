%% Numerical Solution
clear all;
D=1.0;
dt=0.01;
N=128;
m=2;
A=1.0;
kappa=1.0;
nstep=5000;
col_labels={'width_spec';'width_analytical'}; 

% Declarations
conc=zeros(N,1);
c_prime=zeros(N, 1);

% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        conc(i)=1;
    end
    if ((i==N/4) || (i==3*N/4))
        conc(i)=0.5;
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
plot(conc, 'b*');

x1=get_N(conc, 0.1);
x2=get_N(conc, 0.9);
width_spec=x2-x1;


%% Analytical solution

x3=conc_analytical(0.1);
x4=conc_analytical(0.9);

width_analytical=x4-x3;
width_data=[width_spec;width_analytical];

table(col_labels, width_data)

%% Functions

function x_value = conc_analytical(c1)
x_value=2*atanh(2*c1-1);
end

function N_avg = get_N(conc_vector, c_value)
count=0;
sum_N=0;
istep=numel(conc_vector)/2;
for i=1:istep
    if (conc_vector(i)>lim_a && conc_vector(i)<lim_b)
        sum_N=sum_N+i;
        count=count+1;        
    end
end
N_avg=sum_N/count;
end


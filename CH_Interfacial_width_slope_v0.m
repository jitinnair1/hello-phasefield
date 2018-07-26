% In this version, CH solution is by FFT and width calculation is done by taking derivative of conc using central difference followed by intepolation.

clear all;
D=1.0;dt=0.1;N=400;
A=1.0;kappa=1.0;nstep=2000;

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
            k=(i-1)*delk;
        end 
        if ((i) > halfN)
            k=(i-1-N)*delk;
        end
        k2=k*k;
        k4=k2*k2;
        
        c_hat(i)=(c_hat(i)-dt*k2*g_hat(i))/(1+2*k4*dt);
    end
    conc=real(ifft(c_hat));

end
%% Calculating c_prime and interfacial width
c_prime=get_diff(conc, 1.0);

[slope_val, slope_ind]=get_slope_val(c_prime);
width_spec=1/slope_val;
display(width_spec)
%% Functions

% Get c_prime when input is conc
function diff_vector = get_diff(conc_vector, dx)
N=numel(conc_vector);
diff_vector=zeros(1,N);
for i=1:N
    w=i-1;
    e=i+1;
    if(w<1)
        w=w+N;
    end
    if(e>N)
        e=e-N;
    end
    diff_vector(i)=(conc_vector(e)-conc_vector(w))/(2*dx);
end
end

% Get max value of slope and index of max_value until N/2
function [slope_val, slope_ind] = get_slope_val(c_prime_vec)
val_max=-inf;
max_ind=0;
N=numel(c_prime_vec)/2; % because calc is only for first half of profile
for i=1:N
    if (c_prime_vec(i)>val_max)
        val_max=c_prime_vec(i);
        max_ind=i;
    end
end
slope_val=val_max;
slope_ind=max_ind;
end

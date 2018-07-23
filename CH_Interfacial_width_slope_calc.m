%% Numerical Solution
clear all;
D=1.0;dt=0.01;N=128;
A=1.0;kappa=1.0;nstep=8000;
col_labels={'width_spec';'width_analytical'}; 

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
plot(conc,'*b')
title('Composition vs. Distance')

%% Calculating c_prime and interfacial width
c_prime=get_diff(conc, 1);
clf;
plot(c_prime);
ylabel('$\displaystyle\frac{dc}{dx}$','interpreter','latex')
xlabel('N'), title('Variation of slope')
[slope_val, slope_ind]=get_slope_val(c_prime);
y0=conc(slope_ind);
x0=slope_ind;
m=slope_val;
c0=y0-m*x0;
width_spec=get_width(m, c0, 1, 0);

%% Analytical solution
x3=conc_analytical(0.1);
x4=conc_analytical(0.9);
width_analytical=x4-x3;

% Collecting Interfacial results
width_data=[width_spec;width_analytical];
table(col_labels, width_data)

%% Functions

% Solve Analytical Equation
function x_value = conc_analytical(c1)
x_value=2*atanh(2*c1-1);
end

% Get c_prime when input is conc
function diff_vector = get_diff(conc_vector, dx)
N=numel(conc_vector);
diff_vector=zeros(1,N);
for i=1:N
    w=i-1;
    e=i+1;
    if(w<1)
        w=w+1;
    end
    if(e>N)
        e=e-1;
    end
    diff_vector(i)=(conc_vector(e)-conc_vector(w))/2*dx;
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

% calculate interfacial width
function width=get_width(slope, const, lim_max, lim_min)
x1=(lim_min-const)/slope;
x2=(lim_max-const)/slope;
width=x2-x1;
end

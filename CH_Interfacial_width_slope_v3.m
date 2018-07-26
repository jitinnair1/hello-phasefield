%In this version, CH solution is by FDM and width calculation is done by taking
%derivative of conc using central difference followed by intepolation.

clear all;
D=1.0;dx=1;dt=0.1;N=400;A=1.0;kappa=1.0;nstep=2000;
beta1=dt/(dx*dx);beta2=(2*kappa*beta1)/(dx*dx);

% Declarations
conc=zeros(N,1);
conc_old=zeros(N,1);
% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        conc_old(i)=1;
    end
    if ((i==N/4) || (i==3*N/4))
        conc_old(i)=0.5;
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
        conc(i)=conc_old(i) + beta1*(g(w)-2*g(i)+g(e))...
            - beta2*(conc_old(ww)-4*conc_old(w)+6*conc_old(i)...
            -4*conc_old(e) + conc_old(ee));
        conc_old(i)=conc(i);     
    end
end

% plot(conc_old);title('Composition vs. Distance')

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


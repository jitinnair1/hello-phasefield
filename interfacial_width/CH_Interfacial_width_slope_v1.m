%In this version, CH solution is by FFT and width calculation is by taking
%derivative of conc using central difference followed by intepolation at
%ends of the interface. This is a wrong approach. The interpolation should
%be done about the point of symmetry as in v0 and v2 of the code.

clear all;
D=1.0;dt=0.01;N=400;
A=1.0;kappa=1.0;nstep=8000;
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
display(width_spec)

%% Functions

function x_value = get_N(conc_vector, c_value)
istep=numel(conc_vector)/2;
conc_vector=conc_vector-c_value;

for i=1:istep
    if(conc_vector(i)>0)
        x1=i;
        x2=i-1;
        break;
    end
end

N1=conc_vector(x1)+c_value;
N2=conc_vector(x2)+c_value;

m=(N2-N1)/(x2-x1);
c=N1-m*x1;

x_value=(c_value-c)/m;

end


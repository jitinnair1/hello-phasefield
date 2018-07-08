D=1.0;
dx=0.5;
dt=0.001;
N=128;
m=4;
kappa=1.0;
A=1.0;
beta1=dt/dx*dx;
beta2=2*kappa*beta1/dx*dx;

% Declarations
conc=zeros(N,1);
conc_old=zeros(N,1);
c_prime=zeros(N,1);

% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        conc_old(i)=1;
    end
end

plot(conc_old, 'r*')
hold on

% Define g
g=zeros(N,1);
for i=1:N
    g(i)=2*A*conc_old(i).*(1-conc_old(i)).*(1-2*conc_old(i));
end

% Evolve the profile

for j=1:2
    for k=1:10000
        for i=1:N
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
    plot(conc);
end


% Calculation of Interfacial energy
energy1=0;
energy2=0;

for i=1:N
    energy1=energy1 + A*conc(i).*conc(i).*(1-conc(i).*conc(i)).*(1-conc(i).*conc(i));
end 

% Interfacial For FDM?
for i=1:N
    w=i-1;
    e=i+1;
    if (w<1)
        w=N-w;
    end
    if (e>N)
        e=e-N;
    end
    c_prime(i)=(conc(w)-2*conc(i)+conc(e))/dx*dx;
end

for i=1:N
    energy2 = energy2 + kappa*c_prime(i)*c_prime(i);
end

ans1=0.5*energy1;
ans2=0.5*energy2;
ans3=0.5*(energy1 + energy2);

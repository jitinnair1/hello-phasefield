kappa=1.0;
A=1.0;
conc=zeros(N,1);

%% Initial profile
m=2;
N=64;
for i=1:N
    conc(i)=0.5*(1+sin(2*pi*m*i/N));
end

plot(conc, 'r*')
hold on

%% Plot the variation of kappa
for i=0.5:0.5:2.0
        CH_FDM(conc, N, A, i);
        hold on
end
title('Effect of varying kappa');
legend('Initial Profile','kappa=0.5','kappa=1','kappa=1.5','kappa=2');
hold off


%% CH_FDM funtion
function CH_FDM(conc_old, N, A, kappa)
D=1.0;
dx=0.5;
dt=0.01;
m=2;
beta1=dt/dx*dx;
beta2=2*kappa*beta1/dx*dx;

% Declarations
conc=zeros(N,1);
g=zeros(N,1);

% Evolve the profile
for k=1:6000
        
        % Define g
        g=2*A*conc_old.*(1-conc_old).*(1-2*conc_old);
        
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
    xlabel('Distance'), ylabel('Composition')
end

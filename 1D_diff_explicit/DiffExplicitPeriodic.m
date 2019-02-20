% This program solves a sinusoidal diffusion profile using periodic 
% boundary conditions 

clear all;

% Parameters
N=128;
dx=0.5;
dt=0.01;
D=1.0;
m=2; % number of half-waves
alpha=D*dt/((dx)^2);

% Declarartions
conc=zeros(N, 1);
conc_old=zeros(N, 1);

% Initial Conditions
for i=1:N
    
    % Algebraic modification is to bring sin in the range [0, 1] 
    % by adding 1 and multiplying by 0.5 
    
    conc_old(i)=0.5*(1+sin(2*pi*m*i*dx/N));
    
end

plot(conc_old,'b*'); ylabel('Composition'), xlabel('Distance');
title('1D Diffusion Profile')
hold on

for j=1:20
    for k=1:1500
        for i=1:N
            
            % Implementation of Periodic Boundary Conditions
            w=i-1;
            e=i+1;
            if (w==0)
                w=w+N;
            end
            if (e==N+1)
                e=e-N;
            end
            
            % Solve for concenration
            conc(i)=conc_old(i)*(1-2*alpha)+alpha*(conc_old(w)+conc_old(e));
            conc_old(i)=conc(i);
            
        end % end i-loop
    end % end k-loop
    
    plot(conc_old) % Plot every `k_limit` steps
    
end % end j-loop
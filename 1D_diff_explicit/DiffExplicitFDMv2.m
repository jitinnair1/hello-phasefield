% This program plots an the diffusion profile for specified bpundary
% conditions. The `hold on` and `plot` within the loop give stepwise
% evolution of the profile. Use drawnow for ploting an animated version of
% the profile.

clear all;

% Parameters
N=100;
dx=0.5;
dt=0.1;
D=1;
alpha=D*dt/((dx)^2);

% Setting up the system
x=1:N;
conc=zeros(1, N);

% Boundary Condition (Left End)
conc(1)=1; 

% Plot initial profile
plot(x, conc, 'b'), xlabel('Distance'), ylabel('Composition');
title('1D Diffusion Profile');
hold on

% Evolve the profile
for k=1:20
    for j=1:500
        for i=2:N-1
            conc(i)=conc(i)*(1-2*alpha) + alpha*(conc(i-1)+conc(i+1));
        end
        
        % Application of boundary condition for right end included here
        conc(N)=conc(N)*(1-2*alpha) + 2*alpha*(conc(N-1));
    end
    plot (x, conc);
    hold on;
end


clear all;
%Parameters
N=100;
dx=0.5;
dt=0.1;
D=1;
x=1:N;
alpha=D*dt/((dx)^2);

%Declarations
conc=zeros(1, N);

%Boundary Conditions: Left end - Diricilet, Right end - Neumann 
conc(1)=1;
conc(N)=0;

%Plot initial profile
plot(x, conc, 'b'), xlabel('Distance'), ylabel('Composition');
title('1D Diffusion Profile');
hold on

%Evolve 
for k=1:20
    for j=1:500
        for i=2:N-1
            conc(i)=conc(i)*(1-2*alpha) + alpha*(conc(i-1)+conc(i+1));
        end   
    end
    plot (x, conc);
    hold on;
end


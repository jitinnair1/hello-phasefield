%{
Note: This code is incomplete as of date. Matrix inversion does not work
here since matrix is not square. The problem needs to be solved (maybe) by
Sherman-Morrison technique or any other method.
%}


clear all
%Parameters

N=128;
dx=0.5;
dt=0.01;
D=1.0;
alpha=D*dt/((dx)^2);
m=2;

%Declarartions
conc=zeros(N, 1);
conc_old=zeros(N, 1);
A=zeros(N, N+2);

%Initial Conditions

for i=1:N
    % Algebraic modification to bring sin in the range [0, 1]
    % by adding 1 and multiplying by 0.5
    conc_old(i)=0.5*(1+sin(2*pi*m*i*dx/N));
end

plot(conc_old,'b*'); ylabel('Composition'), xlabel('Distance');
title('1D Diffusion Profile')
hold on

%Putting together the matrix

%First Column
for i=2:N
    A(i,1)=0;
end

%Last Column
for i=1:N-1
    A(i,N+2)=0;
end

%Middle terms
for j=1:N
    A(j,j+1)=1+2*alpha;
end

%First Terms
for j=2:N
    A(j,j)= -alpha;
end

%Last Terms
for j=2:N
    A(j-1,j+1)= -alpha;
end

% Boundary Condition Implementation
A(1,1)=-alpha;
A(N,N+2)=-alpha;

for j=1:20
    for k=1:500
        conc=A \ conc_old;
        for i=1:N
            conc_old(i)=conc(i);
        end
    end %end k-loop
    plot(conc_old)
end %end j-loop
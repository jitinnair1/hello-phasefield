%{
For any value of dx<0.4, the solution is unstable. This critical value of dx
is given by stability analysis. For proof, refer Numerical Analysis
by Richard L. Burden, J. Douglas Faires, Annette M. Burden.

The critical condition is that alpha should be less than or equal to 0.5 where
alpha = D*dt/(dx*dx)
%}

iflag=1;

if (iflag==0)
    dx=[0.5; 1; 1.5];
end

if (iflag==1)
    dx=[0.3; 1; 1.5];
end


N=64;
m=4;
wave_length=N/m;
dt=0.1;

np=numel(dx);
istep=numel(dt);

labels=strings(np, numel(dt));
grid_points=zeros(istep, 1);
for i=1:np
    for j=1:istep
        DiffExplicinFunc(N, dx(i), dt(j), m);
        grid_points(i, j)=N/dx(i);
        labels(i,j)="dx="+string(dx(i));
        hold on
    end
end
title('Effect of varying dx');
legend(labels);
hold off
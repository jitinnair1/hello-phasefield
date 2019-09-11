%{
For any value of dx<0.4, the solution is unstable. This critical value of dx
is given by stability analysis. For proof, refer Numerical Analysis
by Richard L. Burden, J. Douglas Faires, Annette M. Burden.

The critical condition is that alpha should be less than or equal to 0.5 where
alpha = D*dt/(dx*dx)
%}

clear all;
iflag=0;

if (iflag==0)
    dx=[0.3; 0.4; 0.5];
    dt=0.1;
end

if (iflag==1)
    dt=[0.4; 0.5; 0.6];
    dx=1;
end

D=1.0;
N=128;
m=1;
wave_length=N/m;

np=numel(dx);
istep=numel(dt);

labels_dx=strings(np, istep);
labels_dt=strings(np, istep);
grid_points=zeros(istep, 1);
alpha_dx=zeros(np, istep);
alpha_dt=zeros(np, istep);

if (iflag==0)
    % Plot for dx variation
    for i=1:np
        for j=1:istep
            DiffExplicinFunc(N, dx(i), dt(j), m);
            grid_points(i, j)=N/dx(i);
            alpha_dx(i, j)=D*dt(j)/(dx(i)*dx(i));
            labels_dx(i,j)="dx="+string(dx(i))+"  "+...
                "alpha="+string(alpha_dx(i,j));
            hold on
        end
    end
    title('Effect of varying dx');
    legend(labels_dx);
    hold off
end

if (iflag==1)
    % Plot for dt variation
    for i=1:np
        for j=1:istep
            DiffExplicinFunc(N, dx(i), dt(j), m);
            grid_points(i, j)=N/dx(i);
            alpha_dt(i, j)=D*dt(j)/(dx(i)*dx(i));
            labels_dt(i,j)="dt="+string(dt(j))+"  "+...
                "alpha="+string(alpha_dt(i,j));
            hold on
        end
    end
    title('Effect of varying dt');
    legend(labels_dt);
    hold off
end
N=64;
m=8;
dt=[0.1; 0.5];
dx=[0.4; 0.8; 1.2];
np=numel(dx);
istep=numel(dt);

labels=strings(np, numel(dt));
dx_by_dt=zeros(istep, 1);
for i=1:np
    for j=1:istep
    CH_FDM_GridPoint(N, m, dx(i), dt(j))
    dx_by_dt(i,j)=dx(i)/dt(j);
    labels(i,j)="dx="+string(dx(i))+" , dt="+string(dt(j));
    hold on
    end
end
title('Effect of varying dx and dt');
legend(labels);
hold off

dx_by_dt(i,j)=dx_by_dt(i,j)/m;

% Format Results before display
col_labels=reshape(labels, [], 1);
gridpoints_per_wavelength=reshape(dx_by_dt, [], 1);
table(col_labels, gridpoints_per_wavelength)
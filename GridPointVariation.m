N=64;
m=2;
dt=0.01;
for i=0.5:0.5:2.0
    CH_FDM_GridPoint(N, m, i, dt)
    hold on
end
title('Effect of varying dx');
legend('dx=0.5','dx=1','dx=1.5','dx=2');
hold off
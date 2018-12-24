%{ 
This alpha is omega/R*T which is non dimensional 
where omega is the interaction parameter
%}

%Using simple for loop
for alpha=1:0.5:3.5
x=0.001:0.001:0.999;
DS=x.*log(x)+(1.-x).*log(1.-x);
DH=alpha.*x.*(1.-x);
DG=DH+DS;
plot(x, DG)
hold on
end
title('Variation of $\Delta G^{\prime} $ with $\frac{\Omega}{RT}$ ','Interpreter','latex', 'fontsize', 24);
xlabel('Composition', 'Interpreter','latex', 'fontsize', 18)
ylabel('$\Delta G^{\prime}$','Interpreter','latex', 'fontsize', 24);
legend('1.0', '1.5', '2.0', '2.5', '3.0', '3.5' )
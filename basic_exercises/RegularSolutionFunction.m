%Using function
x=0.001:0.001:0.999;

for a=0.5:0.5:3.5
plot(x, G(a,x))
hold on
end 

title('$\Delta G$ variation with $\alpha$','Interpreter','latex', 'fontsize', 24);
xlabel('Composition', 'Interpreter','latex', 'fontsize', 18)
ylabel('$\Delta G$','Interpreter','latex', 'fontsize', 24);

function y = G(a, x)
y = 0.15*x + 0.11*(1.-x) + a.*x.*(1.-x) + x.*log(x) + (1.-x).*log(1.-x);
end

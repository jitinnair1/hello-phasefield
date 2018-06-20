%{ 
This alpha is omega/R*T which is non dimensional 
where omega is the interaction parameter
%}

%Using simple for loop
for alpha=0.1:0.9:5
x=0.001:0.001:0.999;
DS=x.*log(x)+(1.-x).*log(1.-x);
DH=alpha.*x.*(1.-x);
DG=DH+DS;
plot(x, DG)
hold on
end
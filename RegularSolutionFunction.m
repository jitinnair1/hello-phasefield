%Using function
x=0.001:0.001:0.999;

for a=0.1:0.5:5
plot(x, G(a,x))
hold on
end 

function y = G(a, x)
y = a.*x.*(1.-x) + x.*log(x) + (1.-x).*log(1.-x);
end

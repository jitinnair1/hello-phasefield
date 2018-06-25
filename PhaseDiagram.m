clear all;

%Declarations
x1=0.01:0.01:0.49;
x2=0.51:0.01:0.99;
alpha=2:0.01:3.5; %for alpha < 2, spinoidal will not happen
N=numel(alpha);
min1=zeros(1,N);
min2=zeros(1,N);

% Find Minima for different alpha values
for i = 1:numel(alpha)
min1(i)=fminbnd(@(x1) G(alpha(i), x1),0.001,0.499, optimset('TolX', 1.e-8));
min2(i)=fminbnd(@(x2) G(alpha(i), x2),0.501,0.999, optimset('TolX', 1.e-8));
end

T=1./alpha;
plot(min1, T),ylabel('Temperature'),xlabel('Composition');
title('Phase Diagram')
hold on
plot(min2, T);

%Free Energy Function
function y=G(alpha, x)
y=alpha.*x.*(1.-x) + x.*log(x)+(1.-x).*log(1.-x);
end







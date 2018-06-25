clear all;

%Declarations
x=0.01:0.01:0.99;
alpha=2.01:0.01:3.5; %for alpha < 2, spinodal will not happen
N=numel(alpha);
min1=zeros(1,N);
min2=zeros(1,N);
root1=zeros(1,N);
root2=zeros(1,N);
x0=[0.001; 0.501];
x1=[0.499; 0.999];

% Find Minima for different alpha values
for i = 1:numel(alpha)
min1(i)=fminbnd(@(x) G(alpha(i), x),0.001,0.499, optimset('TolX', 1.e-8));
min2(i)=fminbnd(@(x) G(alpha(i), x),0.501,0.999, optimset('TolX', 1.e-8));
end

%Find Spinodal points
for i = 1:numel(alpha)
root1(i)=fzero(@(x) Gdprime(alpha(i), x), x0); 
root2(i)=fzero(@(x) Gdprime(alpha(i), x), x1);
end

T=1./alpha;
plot(min1, T),ylabel('Temperature'),xlabel('Composition');
title('Phase Diagram')
hold on
plot(min2, T);
hold on
plot(root1, T);
hold on
plot(root2, T);

%Free Energy Function
function y=G(alpha, x)
y=alpha.*x.*(1.-x) + x.*log(x)+(1.-x).*log(1.-x);
end

%Free Energy Derivative Function
function y=Gdprime(alpha, x)
y=-2*alpha + 1./x +1./(1.-x);
end
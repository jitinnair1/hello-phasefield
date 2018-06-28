% Simple program to check how the cuve varies for different 
% values of A in f = A*c^2(1-c)^2

c=-0.1:0.01:1.1;
A=0.5:0.5:2;
for i=1:numel(A)
    f=A(i).*c.*c.*(1-c).*(1-c);
    plot(c, f); 
    hold on
end
xlabel('Values of c'), ylabel('Value of f');
title('Variation of f with A');
legend('A=0.5','A=1','A=1.5','A=2');
    
% Simple program to check how the cuve varies for different
% values of A in f = A*c^2(1-c)^2

%Physical significance of A is that it gives barrier height

c=-0.1:0.01:1.1;
A=0.5:0.5:2;
for i=1:numel(A)
    f=A(i).*c.*c.*(1-c).*(1-c);
    subplot(2,2,1);
    plot(c, f);
    hold on
end
xlabel('Values of c'), ylabel('Value of f');
title('Variation of f with A');
legend('A=0.5','A=1','A=1.5','A=2');

% Plotting f' (1st derivative)

c=-0.1:0.01:1.1;
A=0.5:0.5:2;
for i=1:numel(A)
    f_prime=2*A(i).*c.*(1-c).*(1-2.*c);
    subplot(2,2,2)
    plot(c, f_prime);
    hold on
end
xlabel('Values of c'), ylabel('Value of f_prime');
title('Variation of f prime with A');
legend('A=0.5','A=1','A=1.5','A=2');


% Plotting f'' (2nd derivative)

c=-0.1:0.01:1.1;
A=0.5:0.5:2;
for i=1:numel(A)
    f_dprime=2*A(i).*(1-c).*(1-2.*c)-2*A(i).*c.*(1-2.*c)-4*A(i).*c.*(1-c);
    subplot(2,2,3)
    plot(c, f_dprime);
    hold on
end
xlabel('Values of c'), ylabel('Value of f double prime');
title('Variation of f double prime with A');
legend('A=0.5','A=1','A=1.5','A=2');

% Combined plot of f, f_prime and f_double prime for A=1

f=A(1).*c.*c.*(1-c).*(1-c);
f_prime=2*A(1).*c.*(1-c).*(1-2.*c);
f_dprime=2*A(1).*(1-c).*(1-2.*c)-2*A(i).*c.*(1-2.*c)-4*A(i).*c.*(1-c);
subplot(2,2,4);
plot(c, f, c, f_prime, c, f_dprime);
xlabel('Values of c'), ylabel('Value of f, f prime and f double prime');
title('Variation of f and derivatives with c');
legend('f','f prime','f double prime');


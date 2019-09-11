%{
error = log(factorial(n)) - n.*log(n) + n
plot error vs. n
%}

n=2:190;
error = (log(factorial(n)) - n.*log(n) + n)./(log(factorial(n)));
plot(n, error), xlabel('Range'), ylabel('Relative Error')
title('Error in Stirling Approximation with N')


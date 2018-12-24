function [dfdeta] = free_energy_diff(i, j, ngrain, etas, eta, igrain)
format long;
A=1.0;
B=1.0;
sum=0.0;

for jgrain=1:ngrain
    if (jgrain ~= igrain)
        sum=sum+etas(i, j, jgrain)*etas(i, j, jgrain);
    end
end

dfdeta=A*(2.0*B*eta(i, j)*sum + eta(i,j)^3 - eta(i, j));

end

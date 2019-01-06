function [dfcon, dfeta] = free_energy_diff(i, j, con, eta, etas, npart, iflag)
format long
A=16.0;
B=1.0;

dfcon=0.0;
dfeta=0.0;

% derivative of free energy for concentration
if (iflag==1)
    sum2=0.0;
    sum3=0.0;
    for ipart=1:npart
        sum2=sum2 + etas(i, j, ipart)^2;
        sum3=sum3 + etas(i, j, ipart)^2;
    end
    
    dfcon=B*(2.0*con(i, j) + 4.*sum3 - 6.0.*sum2) - 2.0*A*con(i, j)^2* ...
        (1.0-con(i, j))+2.0*A*con(i, j)*(1.0-con(i, j))^2;
    
end

% derivative of free energy for etas
if (iflag==2)
    sum2=0.0;
    for ipart=1:npart
        sum2=sum2 + etas(i, j, ipart)^2;
    end
    
    dfeta=B*(-12.0*eta(i, j)^2.*(2.0-con(i, j)) + 12.0*eta(i, j)*...
        (1.0 - con(i, j)) + 12.0*eta(i, j)*sum2);
    
end % if

end % function

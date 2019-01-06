function [npart,etas,con] = micro_sint_pre(Nx,Ny,npart,iflag)

format long;

%--- initialize:

con=zeros(Nx, Ny);
etas=zeros(Nx, Ny, npart); % for all grains
end % if


%---
R1 = 0.2*Nx;
R2 = 0.5*R1;


x1 = Nx/2;
y1 = 0.4*Ny;
y2 = 0.7*Ny;


for i=1:Nx
for j=1:Ny

xx1 =sqrt((i-x1)^2 +(j-y1)^2);
xx2 =sqrt((i-x1)^2 +(j-y2)^2);

if( xx1 <= R1)
con(i,j)=0.999;
etas(i,j,1) =0.999;
end % if

if(xx2 <= R2)
con(i,j) = 0.999;
etas(i,j,1)=0.0;
etas(i,j,2)=0.9999;
end % if

end % j
end % i

end %end function

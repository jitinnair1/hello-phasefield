function [npart,etas,con] = micro_sint_pre(Nx,Ny,npart,sflag)
% initialize:

con=zeros(Nx, Ny);
etas=zeros(Nx, Ny, npart); % for all grains

% sflag==1 for two unequal spherical grains
if(sflag==1)
    
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
    
end

% sflag==2 for two equal hexagonal grains
if (sflag==2)

% Set layout
a=0.2*Nx;
x1 = Nx/2;
y1 = 0.4*Ny;
y2 = 0.7*Ny;
    
% for hexagon 1
x_hex1=a*[-1 -0.5 0.5 1 0.5 -0.5 -1]+x1;
y_hex1=a*sqrt(3)*[0 -0.5 -0.5 0 0.5 0.5 0]+y1;

% for hexagon 2
x_hex2=a*[-1 -0.5 0.5 1 0.5 -0.5 -1]+x1;
y_hex2=a*sqrt(3)*[0 -0.5 -0.5 0 0.5 0.5 0]+y2;


end

% sflag==3 for two equal spherical grains
if (sflag==3)
    
    %---
    R1 = 0.2*Nx;
    R2 = R1;
    
    
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
    
    
end


end %end function

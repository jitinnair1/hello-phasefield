function [etas, ngrain, glist] = initial_struc(Nx, Ny)
% Open Voronoi results
in=fopen('grain_25.inp','r');

% Declarations
twopi=2*pi; epsilon=1.0e-4; ndime=2; %number of dimensions - to be changed for 3D
dx=0.5;dy=0.5;

nvpoin=fscanf(in,'%d', 1);nvnode=fscanf(in,'%d', 1);
nvelem=fscanf(in,'%d', 1);ngrain=fscanf(in,'%d', 1);

for ipoin=1:nvpoin
  jpoin=fscanf(in,'%d',1);
  dummy=fscanf(in, '%lf %lf',[2, 1]);
  %vcord=zeros(jpoin, ndime);
  for idime=1:ndime
    vcord(jpoin, idime) = dummy(idime);
  end
end

for ielem=1:nvelem
  jelem=fscanf(in, '%d', 1);
  dummy=fscanf(in, '%d', [nvnode+1, 1] );
  %vlnods=zeros(ielem,nvnode+1);
  for inode=1:nvnode+1
    vlnods(ielem, inode)=dummy(inode);
  end
end

%nnode2=zeros(nvelem);
for ielem=1:nvelem
  jnode=0;
  for inode=1:nvnode
    knode=vlnods(ielem, inode);
    if(knode~=0)
      jnode=jnode+1;
    end
  end
  nnode2(ielem)=jnode;
end

%gx=zeros(Nx);gy=zeros(Ny);

% Form the Grid
for i=1:Nx
  gx(i)=i*dx;
end

for j=1:Ny
  gy(j)=j*dy;
end

%etas=zeros(Nx, Ny, ngrain);
% Initialise order parameters
for i=1:Nx
  for j=1:Ny
    for igrain=1:ngrain
      etas(i, j, igrain)=0.0;
    end
  end
end

for i=1:Nx
  for j=1:Ny
    for ielem=1:nvelem
      igrain=vlnods(ielem, nvnode+1);
      theta=0.0;
      mnode=nnode2(ielem);
      for inode=1:mnode
        knode=vlnods(ielem, inode);

        xv1=vcord(knode, 1);yv1=vcord(knode, 2);
        jnode=vlnods(ielem, inode+1);

        if (inode==mnode)
          jnode=vlnods(ielem, 1);
        end

        xv2=vcord(jnode, 1);    yv2=vcord(jnode, 2);
        p1x=xv1-gx(i);  p1y=yv1-gy(j);
        p2x=xv2-gx(i);  p2y=yv2-gy(j);
        x1=sqrt(p1x*p1x+p1y*p1y);   x2=sqrt(p2x*p2x+p2y*p2y);

        if (x1*x2<= epsilon)
          theta=twopi;
        else
          tx1=((p1x*p2x+p1y*p2y)/(x1*x2));
        end

        if (abs(tx1)>= 1.0)
          tx1=0.9999999999;
        end

        theta=theta+acos(tx1);
      end
      if (abs(theta-twopi) <= epsilon)
        etas(i, j, igrain)=1.0;
      end
    end
  end
end

% Initialise glist
%glist=zeros(1, ngrain);
for igrain=1:ngrain
    glist(igrain)=1.0;
end
end

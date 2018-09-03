% Set up files
out=fopen('Voroni_vertices.out','w');
out1=fopen('plot_1.out','w');
out2=fopen('final_plot.p','w');
out3=fopen('original_points.out','w');
out4=fopen('cell_1.out','w');
out5=fopen('grain_25.inp','w');

% Parameters
npoin=25;
xmax=32.0;
ymax=32.0;
x0=0.0;
y0=0.0;
extra=2.0;

% Generate Random Points
rng('default');
rng(7);

x=xmax*rand(npoin, 1);
y=ymax*rand(npoin, 1);

% Duplicate the points for creating nine sectors

for ipoin=1:npoin
    jpoin=npoin+ipoin;
    x(jpoin)=x(ipoin);
    y(jpoin)=y(ipoin)-ymax;
end

for ipoin=1:npoin
    jpoin=npoin*2+ipoin;
    x(jpoin)=x(ipoin)+xmax;
    y(jpoin)=y(ipoin)-ymax;
end

for ipoin=1:npoin
    jpoin=npoin*3+ipoin;
    x(jpoin)=x(ipoin)+xmax;
    y(jpoin)=y(ipoin);
end

for ipoin=1:npoin
    jpoin=npoin*4+ipoin;
    x(jpoin)=x(ipoin)+xmax;
    y(jpoin)=y(ipoin)+ymax;
end

for ipoin=1:npoin
    jpoin=npoin*5+ipoin;
    x(jpoin)=x(ipoin);
    y(jpoin)=y(ipoin)+ymax;
end

for ipoin=1:npoin
    jpoin=npoin*6+ipoin;
    x(jpoin)=x(ipoin)-xmax;
    y(jpoin)=y(ipoin)+ymax;
end

for ipoin=1:npoin
    jpoin=npoin*7+ipoin;
    x(jpoin)=x(ipoin)-xmax;
    y(jpoin)=y(ipoin);
end

for ipoin=1:npoin
    jpoin=npoin*8+ipoin;
    x(jpoin)=x(ipoin)-xmax;
    y(jpoin)=y(ipoin)-ymax;
end

% Print original points
for i=1:ipoin
    fprintf(out3, '%4.6e %14.6e \n', x(i), y(i));
end

% Print simulation cell
fprintf(out4,'%14.6e %14.6e\n',x0,y0);	
fprintf(out4,'%14.6e %14.6e\n',xmax,y0);
fprintf(out4,'%14.6e %14.6e\n',xmax,ymax);		
fprintf(out4,'%14.6e %14.6e\n',x0,ymax);		
fprintf(out4,'%14.6e %14.6e\n',x0,y0);

% Generate Voronoi diagram
[c, f] = voronoin([x, y]);

% Rearrange output
nvelem=size(f);
ncount=0;
for i=1:nvelem
    flag=1;
    vnodes=f{i,:};
    nnode=size(vnodes, 2);
    
    for j=1:nnode
        if(vnodes(j)==1)
            flag=0;
        end
    end
    
    if(flag==1)
        ncount=ncount+1;
        for j=1:nnode
            lnods(ncount, j)=vnodes(j);
        end
    end
end

% Print Voronoi results
for i=1:ncount
    fprintf(out, '# i %d\n', i);
    nnode=size(lnods, 2);
    
    for j=1:nnode
        kk=lnods(i, j);
        if(kk ~= 0)
            fprintf(out,'%14.6e %14.6e\n',c(kk,1),c(kk,2));
        end
    end
    kk=lnods(i, 1);
    fprintf(out, '%14.6e %14.6e\n', c(kk, 1), c(kk, 2));
    fprintf(out, '\n');
end

% Clip simulation cell
nelem=0;
for i=1:ncount
    flag=0;
    for j=1:nnode
        kk=lnods(i, j);
        if(kk ~= 0)
            if(c(kk, 1) >= -extra && c(kk, 1) <= xmax+extra)
                if (c(kk, 2) >= -extra && c(kk, 2) <= ymax+extra)
                    flag=1;
                end
            end
        end
    end
    
    if(flag==1)
        nelem=nelem+1;
        jnode=0;
        for j=1:nnode
            kk=lnods(i, j);
            if(kk ~= 0)
                jnode=jnode+1;
                lnods2(nelem, jnode)=lnods(i, j);
            end
        end
        nnode2(nelem)=jnode;
    end
end

% Assign Grain numbers
twopi=2*pi;
epsilon=1.0e-4;
igrain=zeros(nelem, 1);

for isector=1:9
    for ipoin=1:npoin
        jpoin=(isector-1)*npoin+ipoin;
        for ielem=1:nelem
            theta=0.0;
            nnode=nnode2(ielem);
            
            for inode=1:nnode
                kk=lnods2(ielem, inode);
                xv1=c(kk, 1);
                yv1=c(kk, 2);
                
                jnode=inode+1;
                if(inode==nnode)
                    jnode=1;
                end
                
                jj=lnods2(ielem, jnode);
                xv2=c(jj, 1);
                yv2=c(jj, 2);
                
                p2x=(xv1-x(jpoin));
                p2y=(yv1-y(jpoin));
                
                p1x=(xv2-x(jpoin));
                p1y=(yv2-y(jpoin));
                
                x1=sqrt(p1x*p1x+p1y*p1y);
                x2=sqrt(p2x*p2x+p2y*p2y);
                
                if(x1*x2<= epsilon)
                    theta=twopi;
                else
                    tx1=((p1x*p2x+p1y*p2y)/(x1*x2));
                    
                    if (abs(tx1) >= 1.0)
                        tx1=0.9999999999;
                    end
                    theta=theta+acos(tx1);
                end
            end
            
            if( abs(theta-twopi) <= epsilon)
                igrain(ielem)=ipoin;
            end
            
        end
    end
end

% Print input files
[nn1, nn2] = size(c);
nnode=size(lnods2, 2);
fprintf(out5, '%5d %5d %5d %5d\n', (nn1-1), nnode, nelem, npoin);

for i=2:nn1
    fprintf(out5, '%5d %14.6e %14.6e\n', i, c(i, 1), c(i, 2));
end

for i=1:nelem
    fprintf(out5, '%5d', i);
    
    for j=1:nnode
        fprintf(out5, '%5d', lnods2(i, j));
    end
    
    fprintf(out5,'%5d', igrain(i));
    fprintf(out5, '\n');
end

for i=1:nelem
	fprintf(out1,'# i %d %d\n',i,nnode2(i));

nnode=size(lnods2,2);

ncount=0;
xcod =0.0;
ycod =0.0;

%nnode = nnode2(ielem);

for j=1:nnode

kk = lnods2(i,j);
fprintf(out1,'# %5d %5d\n',j,kk);

if(kk ~= 0)
fprintf(out1,'%14.6e %14.6e\n',c(kk,1),c(kk,2));
ncount=ncount+1;
xcod =xcod +c(kk,1);
ycod =ycod +c(kk,2);
end
end
kk = lnods2(i,1);
fprintf(out1,'%14.6e %14.6e\n',c(kk,1),c(kk,2));
fprintf(out1,'\n');

xcod =xcod/ncount;
ycod =ycod/ncount;

fprintf(out2,'set label  ');fprintf(out2,'"'); fprintf(out2,'%d',i); ...
fprintf(out2,'" at'); fprintf(out2,'%14.6e  , %14.6e\n',xcod,ycod);
fprintf(out2,'\n');

fprintf(out2,'set label  ');fprintf(out2,'"'); fprintf(out2,'%d',igrain(i)); ...
fprintf(out2,'" at'); fprintf(out2,'%14.6e  , %14.6e\n',xcod+1.5,ycod+1.5);
fprintf(out2,'\n');

end
fprintf(out2, 'plot "plot_1.out" w l, "cell_1.out" w l, "original_points.out"\n');
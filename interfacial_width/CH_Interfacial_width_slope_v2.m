%In this version, CH solution is by FDM and width calculation is done 
%by taking simple intepolation.

clear all;
D=1.0;dx=1;dt=0.1;N=800;A=1.0;kappa=1.0;nstep=2000;
beta1=dt/(dx*dx);beta2=(2*kappa*beta1)/(dx*dx);

% Declarations
conc=zeros(N,1);
conc_old=zeros(N,1);
% Initial profile
for i=1:N
    if (i>N/4 && i<3*N/4)
        conc_old(i)=1;
    end
    if ((i==N/4) || (i==3*N/4))
        conc_old(i)=0.5;
    end
end
% Define g
g=zeros(N,1);
% Evolve the profile
for k=1:nstep
    for i=1:N
        % Define g
        g=2*A*conc_old.*(1-conc_old).*(1-2*conc_old);
        w=i-1;
        ww=i-2;
        e=i+1;
        ee=i+2;
        if (ww<1)
            ww=ww+N;
        end
        if (w<1)
            w=w+N;
        end
        if (ee>N)
            ee=ee-N;
        end
        if(e>N)
            e=e-N;
        end
        % ellipsis are for text wrap in MATLAB editor window
        conc(i)=conc_old(i) + beta1*(g(w)-2*g(i)+g(e))...
            - beta2*(conc_old(ww)-4*conc_old(w)+6*conc_old(i)...
            -4*conc_old(e) + conc_old(ee));
        conc_old(i)=conc(i);     
    end
end

% plot(conc_old);title('Composition vs. Distance')

%% Calculating slope
lim=N/4;lima=lim-1;limb=lim+1;
m=(conc_old(limb)-conc_old(lima))/(limb-lima);
width_spec=1/m;

% Collecting Interfacial results
error_width=((width_spec-4.0)/4.0);
display(width_spec);display(error_width);

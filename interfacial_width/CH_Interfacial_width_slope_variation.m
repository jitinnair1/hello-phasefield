% Code plots variation of IF width on changing number of timesteps i.e.
% "nsteps". Refer CH_Interfacial_width_func.m for details of the function used.

%In the function called here, the solution of CH is done by FFT and calculation of slope is
%done by simple linear interpolation about the point of symmetry

clear all;
nstep=200:200:6000;
count=numel(nstep);
width=zeros(1, count);

for i=1:count
    width(i)=CH_Interfacial_width_func(nstep(i));
end

plot(nstep, width);
title('Variation of width with number of time steps')
ylabel('Width'), xlabel('No. of steps')

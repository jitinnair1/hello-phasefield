clear all;
N=128;
D=1;
m=2;
dt=0.1;

%Declararions
conc=zeros(N,1);

%Initial Profile
for i=1:N
    conc(i)=sin(2*pi*m*i/N);
    % conc(i)=0.5*(1+sin(2*pi*m*i/N));
end

%Plot initial profile
plot(conc, 'r*');
hold on

%Periodic Boundary Definitions
halfN=N/2;
delk=2*pi/N;

%Evolution
for m=1:15 %Loop for time steps
    
    conc_hat=fft(conc);
    
    for n=1:300 %Loop for calculation
        
        for i=1:N
            
            if (i<halfN)
                k=i*delk;
            end
            
            if (i>=halfN)
                k=(i-N)*delk;
            end
            
            %Implict condition
            conc_hat(i,1)=conc_hat(i,1)/(1+D*k*k*dt); 
        end
    end
    
    conc=real(ifft(conc_hat));
    plot(conc);
end


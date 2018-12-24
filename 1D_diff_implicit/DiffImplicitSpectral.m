clear all;
N=128;
D=1;
m=2;
dt=0.01;

%Declararions
conc=zeros(N,1);

%Initial Profile
for i=1:N
    %conc(i)=sin(2*pi*m*i/N);
    conc(i)=0.5*(1+sin(2*pi*m*i/N));
end

%Plot initial profile
plot(conc, 'r*');
hold on

%Periodic Boundary Definitions
halfN=N/2;
delk=2*pi/N;

%Evolution
for m=1:8 %Loop for time steps
   
    
    for n=1:1500 %Loop for calculation
        
        
        
        for i=1:N
            
            conc_hat=fft(conc);
            
            % i-1 so as to include the k=0 point
            if ((i-1)<halfN)
                k=(i-1)*delk;
            end
            
            if ((i-1)>=halfN)
                k=(i-1-N)*delk;
            end
            
            %Implict condition
            conc_hat(i,1)=conc_hat(i,1)/(1+D*k*k*dt); 
            
            conc=real(ifft(conc_hat));
        end
    end
    
    
    plot(conc);
end


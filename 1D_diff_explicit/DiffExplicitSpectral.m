clear all;
N=128;
D=1;
m=2;
dt=0.01;
dx=1;
alpha=D*dt/(dx*dx);

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
        
        
        
        for i=1:N-1
            
            conc_hat=fft(conc);
            
            w=i-1;
            e=i+1;
            
            % i-1 so as to include the k=0 point
            if (w<1)
                w=i+N;
            end
            
            if (e>N)
                e=i-N;
            end
            
            %Implict condition
            conc_hat(i,1)=conc_hat(i,1) + alpha*(conc_hat(w,1)...
                - 2*conc_hat(i,1) + conc_hat(e,1)); 
            
            conc=real(ifft(conc_hat));
        end
    end
    
    
    plot(conc);
end


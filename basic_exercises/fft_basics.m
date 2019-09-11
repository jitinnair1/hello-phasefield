% Clear
clear all;
clc
clf

% Number of grid points
N = 16;

% Make an initial sinusoidal profile of composition with M waves
c = zeros(N,1);
M=1;
for j=1:N
c(j) = sin(2*pi*M*j/N);
end

%Plot original concentration profile
subplot(2,2,1)
plot(c)
title('Original Composition')
hold on

%Plot real part
chat = fft(c);
subplot(2,2,2)
plot(real(chat),'r+-')
title('Real Part')

%Plot imaginary part
subplot(2,2,3)
plot(imag(chat),'g+-')
title('Imaginary Part')

% Shift to centre
a = fftshift(chat);
subplot(2,2,4)
plot(imag(a));
title('After Shift')

%% For cos wave

figure()
c = zeros(N,1);
M=1;
for j=1:N
c(j) = cos(2*pi*M*j/N);
end

%Plot original concentration profile
subplot(2,2,1)
title('Original Composition')
plot(c)
hold on

%Plot real part
chat = fft(c);
subplot(2,2,2)
plot(real(chat),'r+-')
title('Real Part')

%Plot imaginary part
subplot(2,2,3)
plot(imag(chat),'g+-')

% Shift to centre
a = fftshift(chat);
subplot(2,2,4)
plot(imag(a));

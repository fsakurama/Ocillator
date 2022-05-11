% Van der Pol Oscillator Force

clear;

fs = 44100; % sample rate (Hz)
T = 1/fs; % time period (s)
dur = 1; % duration of signal (s)
Ns = floor(dur*fs); % number of samples
t = (0:Ns-1)*T; % time vector

beta = 50;
k = 1e4;
m = 0.01;

F = 1e4*triang(Ns);

um = 0; % x(n-1)
u = 0; % x(n)
up = 0; % x(n+1)

out = zeros(1,Ns); % output vector
out2 = zeros(1,Ns);

for n=1:Ns
    
    A = 2*u*(4-(T^2)*k/m) - um*(4+2*T*beta/m+(T^2)*k/m) + 2*T*beta*(u^2)*um/m + 4*(T^2)*F(n)/m;
    B = 4 - 2*T*beta*(1-(u^2))/m + (T^2)*k/m;
    up = A/B;
    
    out(n) = u;
    out2(n) = (up-um)/(2*T);
    
    um = u;
    u = up;
    
end

soundsc(out,fs);

% displacement-time plot
figure(1);
clf;
plot(t,out);
xlabel('time (s)');
ylabel('displacemant (m)');
title('Displacement-Time Plot');
grid;

% Periodic change
figure(2);
clf;
plot(out,out2);
xlabel('u');
ylabel('(up-um)/(2*T)');
title('Periodic change');
grid;
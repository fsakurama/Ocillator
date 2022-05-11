% Linear Oscillator

clear;

fs = 44100; % sample rate (Hz)
T = 1/fs; % time period (s)
dur = 0.01; % duration of signal (s)
Ns = floor(dur*fs); % number of samples
t = (0:Ns-1)*T; % time vector

m = 0.001; % mass (kg), Increasing mass: low frequency
k = 1e4; % spring constant, High spring constant: high frequency

um = 0.01; % u(n-1)
u = 0; % u(n)
up = 0; % u(n+1)

out = zeros(1,Ns); % output vector
KE = zeros(1,Ns); % kinetic energy
PE = zeros(1,Ns); % potential energy

for n=1:(Ns)
    
    up = (2-k*T^2/m) * u - um;
    
    KE(n) = 0.5*m*(((up-u)/T)^2);
    PE(n) = 0.125*k*(up+u)^2;
    
    out(n) = u;
    
    um = u;
    u = up;   
    
end

soundsc(out,fs);

TE = PE + KE; % calculate total energy

% displacement-time plot
figure(1);
clf;
plot(t,out);
xlabel('time (s)');
ylabel('displacemant (m)');
title('Displacement-Time Plot');
grid;

% energy terms plot
figure(2);
clf;
plot(t,TE,'k-','Linewidth',1.5);
hold on;
plot(t,KE,'r-','Linewidth',1.0);
plot(t,PE,'b-','Linewidth',1.0);
hold off;
xlabel('time (s)');
ylabel('energy terms');
title('Energy');
grid;
legend('Total Energy','Kinetic Energy','Potential Energy','Location','northeast');
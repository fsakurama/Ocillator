% Nonlinear Ocscillator

clear;

fs = 44100; % sample rate (Hz)
T = 1/fs; % time period (s)
dur = 1; % duration of signal (s)
Ns = floor(dur*fs); % number of samples
t = (0:Ns-1)*T; % time vector

m = 0.5; % mass (kg)
k3 = 5e3; % spring constant

um = 0.01; % u(n-1)
u = 0; % u(n)
up = 0; % u(n+1)

psim = 0; % psi(n-1)
psi = 0; % psi(n)
psip = 0; % psi(n+1)

g = 0;

out = zeros(1,Ns); % output vector
out2 = zeros(1,Ns); % output vector
KE = zeros(1,Ns); % kinetic energy
PE = zeros(1,Ns); % potential energy

for n=1:(Ns)
    
    g = sqrt(2*k3)*u; 
    up = (m*(2*u-um)/T^2 + g*((g*um/2)-psim)) / (m/(T^2) + (g^2)/2);
    psip = psim + g*(up-um);
    
    KE(n) = 0.5*m*(((up-u)/T)^2);
    PE(n) = 0.125*(psip+psi)^2;
    
    out(n) = u;
    out2(n) = (up-um)/(2*T);
    
    um = u;
    u = up;   
    psim = psi;
    psi = psip;
    
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

% Periodic change
figure(3);
clf;
plot(out,out2);
xlabel('u');
ylabel('(up-um)/(2*T)');
title('Periodic change');
grid;
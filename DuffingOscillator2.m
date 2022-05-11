% Duffing Oscillator 2

clear;

fs = 44100; % sample rate (Hz)
T = 1/fs; % time period (s)
dur = 1; % duration of signal (s)
Ns = floor(dur*fs); % number of samples
t = (0:Ns-1)*T; % time vector

m = 1; % mass (kg)
k1 = 1e4; % linear stiffness
k3 = 1e3; % non-liniarity in the restoring force
r = 3; % damping coefficient, r>0 yields chaos
gamma = 100; % amplitude of periodic driving force
omega = 2*pi*30; % angular frequency
F = gamma*cos(omega*t);

omega0 = sqrt(k1/m); % natural frequency (rad/sec)
alpha = r/(2*m); % decay rate
omegar = sqrt((omega0)^2 - (alpha)^2); % resonance frequency
R = exp(-alpha*T);
xi = (T^2)/m;
a = (1 - 2*R*cos(omegar*T) + (R^2))/(1+ 2*R*cos(omegar*T) + (R^2));
b = (2*(1-(R^2)))/(1 + 2*R*cos(omegar*T) + (R^2));

um = 0.1; % u(n-1)
u = 0; % u(n)
up = 0; % u(n+1)

psim = 0; % psi(n-1)
psi = 0; % psi(n)
psip = 0; % psi(n+1)

g = 0;

out = zeros(1,Ns); % output vector
out2 = zeros(1,Ns); % output vector
out3 = zeros(1,Ns);
out4 = zeros(1,Ns);

KE = zeros(1,Ns); % kinetic energy
PE = zeros(1,Ns); % potential energy

for n=1:(Ns)
    
    g = sqrt(2*k3)*u;
    up = (2*(1-a)*u - (1+a-b)*um + xi*(F(n)+(0.25*g*um-psim)*g)) / (1 + a + 0.25*(g^2)*xi + b);
    psip = psim + g*(up-um);
    
    KE(n) = 0.5*m*(((up-u)/T)^2);
    PE(n) = 0.125*k1*(up+u)^2+0.125*(psip+psi)^2;
    
    out(n) = u;
    out2(n) = (up-um)/(2*T);
    out3(n) = -k1*u;
    out4(n) = -k1*u - k3*(u^3);
    
    um = u;
    u = up;   
    psim = psi;
    psi = psip;
    
end

TE = PE + KE; % calculate total energy

soundsc(out,fs);

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

% restoring force
figure(4);
clf;
plot(out,out4);
hold on;
plot(out,out3);
hold off;
xlabel('u');
ylabel('F');
title('restoring force');
grid;
legend('-k1*u','-k1*u - k3*(u^3)','Location','northeast');
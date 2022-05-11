% Collision 2

clear;

fs = 44100; % sample rate (Hz)
T = 1/fs; % time period (s)
dur = 2; % duration of signal (s)
Ns = floor(dur*fs); % number of samples
t = (0:Ns-1)*T; % time vector

m = 0.01;
kc = 1e5;
alp = 1.5;

um = 0.1; % u(n-1)
u = 0; % u(n)
up = 0; % u(n+1)

up0 = 0;

psim = 0; % psi(n-1)
psi = 0; % psi(n)
psip = 0; % psi(n+1)

k = 3000;
z = 0;
a = k*T^2/m;
xi = (T^2)/m;

out = zeros(1,Ns);
out2 = zeros(1,Ns);
out3 = zeros(1,Ns);

KE = zeros(1,Ns); % kinetic energy
PE = zeros(1,Ns); % potential energy

for n=1:Ns
    
    h = (2-a)*u -2*um;
    up0 = h/(1+a)+um;
    g = 2*(-psim)/((up0-z)-(um-z));
    
    if u-z>0 && kc>0
        g = sqrt(0.5*kc*(alp+1)*(u-z)^(alp-1));
    end 
    
    s = (h-xi*(psim-0.25*g*z)*g) / (1+a+0.25*xi*g^2);
    up = s+um;
    psip = psim+0.5*g*s;
    
    out(n) = u;
    out2(n) = (up-um)/(2*T);
    out3(n) = g;
    
    KE(n) = 0.5*m*(((up-u)/T)^2);
    PE(n) = 0.125*k*(up+u)^2;

    psim = psip;
    um = u;
    u = up;
end

TE = PE + KE;

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

figure(3);
clf;
plot(t,out3);
xlabel('time (s)');
ylabel('g');
title('g');

% energy terms plot
figure(4);
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
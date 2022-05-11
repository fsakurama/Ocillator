% Two-mass Oscillator

clear;

soundtype = 2; % 1:displacement, 2:velocity, 3:force

fs = 44100; % sample rate (Hz)
T = 1/fs; % time period (s)
dur = 2; % duration of signal (s)
Ns = floor(dur*fs); % number of samples
t = (0:Ns-1)*T; % time vector

m1 = 1e-2; % mass (kg)
m2 = 1e-2;
k1 = 1e4; % spring constant
k2 = 1e3;

a1 = 0.25*k1*(T^2)/m1;
a2 = 0.25*k1*(T^2)/m2;
a3 = 0.25*k2*(T^2)/m2;

M0 = [1+a1 -a1; -a2 1+a2+a3];
M1 = [2*(1-a1) 2*a1; 2*a2 2*(1-a2-a3)];
M2 = [-2*(1+a1) 2*a1; 2*a2 -2*(1+a2+a3)];

um = [0.01; 0];
u = [0; 0];

out1 = zeros(1,Ns); % output vector
out2 = zeros(1,Ns);
outg = zeros(1,Ns);
out1v = zeros(1,Ns); % output vector
out2v = zeros(1,Ns); % output vector
outF1 = zeros(1,Ns);
outF2 = zeros(1,Ns);

KE = zeros(1,Ns); % kinetic energy
PE = zeros(1,Ns); % potential energy

for n=1:Ns
    
    s = M0 \ (M1*u + M2*um);
    
    up = s + um;
    
    KE(n) = 0.5*m1*(((up(1)-u(1))/T)^2) + 0.5*m2*(((up(2)-u(2))/T)^2);
    PE(n) = k1*((up(1)-up(2)+u(1)-u(2))^2)/8 + k2*((up(2)+u(2))^2)/8;
    
    out1(n) = u(1);
    out2(n) = u(2);
    out1v(n) = 0.5*(up(1)-um(1))/T;
    out2v(n) = 0.5*(up(2)-um(2))/T;
    outF1(n) = k1*(u(1)-u(2));
    outF2(n) = k2*u(2);
    
    um = u;
    u = up;
    
end

if soundtype == 1 
    y1 = out1;
    y2 = out2;
elseif soundtype == 2
    y1 = out1v;
    y2 = out2v;
else
    y1 = outF1;
    y2 = outF2;
end

soundsc([y1; y2],fs);

TE = PE + KE;

% energy terms plot
figure(1);
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

%figure(2);
%clf;

%subplot(3,1,1);
%plot(t,out1,'b-',t,out2,'r-');
%xlabel('time (s)');
%ylabel('displacement (m)');
%grid;
%legend('m_1','m_2');

%subplot(3,1,2);
%plot(t,outF1,'b-',t,outF2,'r-');
%hold off;
%xlabel('time (s)');
%ylabel('Force (N)');
%grid;
%legend('F_1','F_2');

%subplot(3,1,3);
%plot(t,out1v,'b-',t,out2v,'r-');
%hold off;
%xlabel('time (s)');
%ylabel('velcoity (m/s)');
%grid;
%legend('v_1','v_2');

figure(2);
clf;
plot(t,out1,'b-',t,out2,'r-');
xlabel('time (s)');
ylabel('displacement (m)');
grid;
legend('m_1','m_2');
title('Displacement');

figure(3);
clf;
plot(t,outF1,'b-',t,outF2,'r-');
hold off;
xlabel('time (s)');
ylabel('Force (N)');
grid;
legend('F_1','F_2');
title('Force');

figure(4);
clf;
plot(t,out1v,'b-',t,out2v,'r-');
hold off;
xlabel('time (s)');
ylabel('velcoity (m/s)');
grid;
legend('v_1','v_2');
title('Velocity');
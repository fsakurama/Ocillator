% 2 - mass Collision Oscillator

clear;

soundtype = 2; % 1:displacement, 2:velocity, 3:force

fs = 44100; % sample rate (Hz)
T = 1/fs; % time period (s)
dur = 5; % duration of signal (s)
Ns = floor(dur*fs); % number of samples
t = (0:Ns-1)*T; % time vector

FAC = 1;
k1 = 1e4*FAC;
m1 = 0.013;
k2 = 9.1e3*FAC;
m2 = 0.011;

kc = 1e7;
alp = 1.3;

a1 = 0.25*k1*(T^2)/m1;
a2 = 0.5*(T^2)/m1;
a3 = 0.25*k1*(T^2)/m2;
a4 = 0.25*k2*(T^2)/m2;

EPSIL = 1e-10;

um = [0; 0];
u = [0; 0];
psim = [0; 0];
psi = [0; 0];
psip = [0; 0];

omega = 2*pi*2;
%ub_sig = zeros(1,Ns);
ub_sig = -0.1*cos(omega*t)+0.099;
%ub_sig = -0.1*triang(Ns);
%ub_sig = -0.1*ones(1,Ns)+0.099;
%ub_sig = -0.1*(1-hann(Ns)')+0.09;

ubm = 0;
ub = 0;
ubp = 0;

out1 = zeros(1,Ns); % output vector
out2 = zeros(1,Ns);
outg = zeros(1,Ns);
out1v = zeros(1,Ns);
out2v = zeros(1,Ns);
outF1 = zeros(1,Ns);
outF2 = zeros(1,Ns);

KE = zeros(1,Ns); % kinetic energy
PE = zeros(1,Ns); % potential energy

for n=1:Ns
    
    ub = 1*ub_sig(n);
    
    M00 = [1+a1 -a1; -a3 1+a3+a4];
    M1 = [2*(1-a1) 2*a1; 2*a3 2*(1-a3-a4)];
    M2 = [-2*(1+a1) 2*a1; 2*a3 -2*(1+a3+a4)];
        
    M00inv = (1/((1+a1)*(1+a3+a4)-a1*a3)) * [1+a3+a4 a1; a3 1+a1];
    s0 = M00inv*M1*u + M00inv*M2*um;
        
    if u(1)-ub>0 && kc>0 
        g1 = sqrt(0.5*kc*(alp+1)*(u(1)-ub)^(alp-1));
    else
        g1 = 2*(-psim(1)*s0(1))/(s0(1)^2 + EPSIL);
    end 
        
    g = [g1; 0];
    
    M0 = [1+a1+0.5*(g(1)^2)*a2 -a1; -a3 1+a3+a4];
    M3 = [-a2*(2*psim(1)-0.5*g1*(ubp-ubm)) 0; 0 0];
    
    M0inv = (1/((1+a1+0.5*(g(1)^2)*a2)*(1+a3+a4)-a1*a3)) * [1+a3+a4 a1; a3 1+a1+0.5*(g(1)^2)*a2];
    s = M0inv*M1*u + M0inv*M2*um + M0inv*M3*g;
        
    up = s + um; 
    psip(1) = psim(1) + 0.5*g(1)*(s(1)-ubp+ubm);

    out1(n) = u(1);
    out2(n) = u(2);
    outg(n) = g(1);
    out1v(n) = 0.5*(up(1)-um(1))/T;
    out2v(n) = 0.5*(up(2)-um(2))/T;
    outF1(n) = k1*(u(1)-u(2));
    outF2(n) = k2*u(2);
        
    KE(n) = 0.5*m1*(((up(1)-u(1))/T)^2) + 0.5*m2*(((up(2)-u(2))/T)^2);
    PE(n) = k1*((up(1)-up(2)+u(1)-u(2))^2)/8 + k2*((up(2)+u(2))^2)/8 + 0.5*(psip(1)^2);
    
    um = u;
    u = up;
    psim = psi;
    psi = psip;
    ubm = ub;
    ub = ubp;
        
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
%plot(t,ub,'k--');
%hold on;
%plot(t,out1,'b-',t,out2,'r-');
%hold off;
%xlabel('time (s)');
%ylabel('displacement (m)');
%grid;
%legend('barrier','m_1','m_2');

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
plot(t,ub_sig,'k--');
hold on;
plot(t,out1,'b-',t,out2,'r-');
hold off;
xlabel('time (s)');
ylabel('displacement (m)');
grid;
legend('barrier','m_1','m_2');
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

figure(5);
clf;
plot(t,outg);
xlabel('time (s)');
ylabel('g');
title('g');
grid;
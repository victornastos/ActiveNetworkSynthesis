%AEM: 9297 ON/EP : NASTOS VIKTOR
%HIGH PASS BUTTERWORTH
a1 = 9;
a2 = 2;
a3 = 9;
a4 = 7;
%a4=7 
m=2;

fp = ( 4 + m ) * 1000;
fs = fp / 2.6;
amin = 24+ ( a3 * ( 6 / 9 ) );
amax = 0.5 + ( a4 / 36 );
ws = 2 * pi * fs;
wp = 2 * pi * fp;
Wp = 1;
Ws = wp / ws;

n  = log10 ( ( 10^( amin / 10 ) - 1 ) / ( 10^( amax / 10 ) - 1 ) ) / ( 2 * log10 (Ws) )
n = ceil(n);
Whp = 1 / ( ( 10^ ( amax / 10 ) - 1 )^( 1 / ( 2 * n ) ) );
w0 = wp/ Whp;

p1= -1;
p2 = -0.809 + 0.587 * 1i;
p3 = -0.809 - 0.587 * 1i;
p4 = -0.309 + 0.951 * 1i;
p5 = -0.309 - 0.951 * 1i;
Q1=0.5;
Q23=0.618;
Q45= 1.618;
W1 = 1;
W23=1;
W45=1;

%MONADA I
k1 = 1;
C11=1;
kf=w0;
km= C11/( kf * 10^(-8));
k1=1;
R11=1 * km;

%MONADA II
C21old=1;
C22old=1;
R21old=1;
R22old=1;
k2=3-(1/Q23);
r21=1;
r22=2-(1/Q23);
%KLIMAKOPOIISI MONADAS II
kf2=w0;
C21new=10^(-8);
km2=km;
R21new = km2 * R21old;
R22new = km2 * R22old;
C21new = C21old / ( kf2 * km2 );
C22new = C22old / ( kf2 * km2 );
r21new= km2*r21;
r22new= km2*r22;

%MONADA III
C31old=1;
C32old=1;
R31old=1;
R32old=1;
k3=3-(1/Q45);
r31=1;
r32=2-(1/Q45);
%KLIMAKOPOIISI MONADAS III
kf3=w0;
C31new=10^(-8);
km3=km;
R31new = km3 * R31old;
R32new = km3 * R32old;
C31new = C31old / ( kf3 * km3 );
C32new = C32old / ( kf3 * km3 );
r31new= km3*r31;
r32new= km3*r32;


ktotal=k1*k2*k3;
aGain = 10^( 0.25 ) / ktotal;
Z1=10*10^3;
Z2=aGain*Z1;

T1 = tf([0 k1*1 0],[0 1 (w0)]);
T2 = tf([k2*1 0 0],[1 ( w0 / Q23 ) (w0)^2]);
T3 = tf([k3*1 0 0],[1 ( w0 / Q45 ) (w0)^2]);

THP = T1*T2*T3;

ltiview ({'bodemag'}, THP)
THP = aGain * THP;

ltiview({'bodemag'}, T1)
ltiview({'bodemag'}, T2)
ltiview({'bodemag'}, T3)
ltiview ({'bodemag'}, THP)
ltiview({'bodemag'}, T1, T2,T3, THP)
plot_transfer_function(THP , [fp fs] );
plot_transfer_function(1/THP , [fp fs] );

%eisagw periodiko shma
f11= (0.5*ws)/(2*pi);
f12=(0.8*ws)/(2*pi);
f13= (1.2*wp) /(2*pi);
f14=(2.4*wp) / (2*pi);
f15= (3.5*wp)/(2*pi) ;
Fs= 200*10^3;
T=0.003;
dt=1/Fs;
t=0:dt:(T);
u1= cos(2*pi*f11*t)+0.6* cos(2*pi* f12*t)+1.5*cos(0.5* f13*t)+0.7*cos(2*pi*f14*t)+0.4*cos(2*pi*f15*t);
figure
plot(u1)
N=T/dt;
figure
lsim(THP,u1,t)
xt=lsim(THP,u1,t);

figure()
plot(t,xt)
n=2^nextpow2(N);
xfourier= fft(xt,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
figure()
plot(f,p1)

nfft=n;
y=fft(u1,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*fs/nfft; 
figure()
plot(f,y_mag)

function plot_transfer_function( tf, frequency_markers )
%PLOT_TRANSFER_FUNCTION Plots bode of a transfer function with markers
%
%   tf                - The transfer function (created using tf)
%   frequency_markers - A matrix of frequencies in Hz
%
%   Example:
%       plot_transfer_function( tf([1000], [1 1000]), [10 1000 10000] );

figure;
x_space = logspace(1,5,5000); % 5000 points between 10^1 and 10^5
x_space = 2 * pi * x_space; % to rad / sec
[mag,~,wout] = bode(tf,x_space);
mag = squeeze(mag);
wout = squeeze(wout);
mag = 20*log10(mag);
wout = wout/2/pi;
semilogx(wout,mag,'-b');
axis([min(wout) max(wout) min(mag)-10 max(mag)+10]);
[num,den] = tfdata(tf);
syms s;
d1 = digits(5);
ltx = latex(vpa(poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s)));
digits(d1);
title(strcat('$',ltx,'$'), 'Interpreter','latex', 'FontSize', 24);
xlabel('Frequency (Hz)', 'FontSize', 18);
ylabel('Magnitude (dB)', 'FontSize', 18);
grid on;
hold on;
[dbMarks,~,frequency_markers] = bode(tf,2 * pi * frequency_markers);
dbMarks = squeeze(dbMarks);
frequency_markers = squeeze(frequency_markers);
dbMarks = 20*log10(dbMarks);
frequency_markers = frequency_markers/2/pi;
Aw = cell(size(frequency_markers, 1) + 1, 1);
Aw{1} = 'Transfer function';
for i = 1 : size(frequency_markers, 1)
    semilogx(frequency_markers(i),dbMarks(i),'o');
    Aw{i+1} = sprintf('Attenuation at %.2f Hz is %.2f dB', ...
        frequency_markers(i), dbMarks(i));
end
legend(Aw,'Location','best','FontSize',12);
set(gca,'FontSize',14);
end











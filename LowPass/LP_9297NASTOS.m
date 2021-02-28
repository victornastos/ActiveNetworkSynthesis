
%AEM: 9297 ON/EP : NASTOS VIKTOR
%LOW PASS CHEBYSHEV
a1=9;
a2=2;
a3=9;
a4=7;

m=1;
fp=(3+m)* 1100;
fs=1.6*fp;
amin= 21.5 + (max(1, a3)-5)*(0.5);
amax = 0.55 +(max(1, a4)-5)/10;

wp= 2 *pi * fp;
ws= 2*pi*fs;
WS= ws/wp;
e= sqrt(10^(amax / 10) - 1);
n  = acosh( sqrt( ( 10^(amin / 10) - 1 ) / ( 10^(amax / 10) - 1) ) ) / acosh(WS);
n_1 = ceil(n);
Whp = cosh( (1 / n_1) * acosh( 1 / e) );
whp = 2 * pi * Whp * fp;

a = asinh(1/e)/n;

ps_k1 = 0;
ps_k2 = 36;
ps_k3 = -36;
ps_k4 = 72;
ps_k5= -72;

s_k1 = -sinh(a)*cosd(ps_k1) + j*cosh(a)*sind(ps_k1);
s_k2 = -sinh(a)*cosd(ps_k2) + j*cosh(a)*sind(ps_k2);
s_k3 = -sinh(a)*cosd(ps_k3) + j*cosh(a)*sind(ps_k3);
s_k4 = -sinh(a)*cosd(ps_k4) + j*cosh(a)*sind(ps_k4);
s_k5 = -sinh(a)*cosd(ps_k5) + j*cosh(a)*sind(ps_k5);

W1 = abs(s_k1);
W23 = abs(s_k2);
W45 = abs(s_k4);

Q_1 = W1 / ( -2 * real(s_k1) );
Q_2 = W23 / ( -2 * real(s_k2) );
Q_3 = W45 / ( -2 * real(s_k4) );

w00= W1 * wp;
w01 = W23 * wp;
w02 = W45 * wp;
%STRATHGIKH (3), 0.01μF, 0dB
%MONADA I

R11old=1;
C11old=1;
C11new=10^(-8);
kf1 = w00;
km1=(C11old)/(kf1*C11new);
R11new=km1*R11old;

%MONADA II 
C21old = 1;
C22old = 1/Q_2;
R21old = 1;
R22old = Q_2;

kf2=w01;
km2 = 10^8*C22old / kf2;
R21new=km2*R21old;
R22new=km2*R22old;
C22new=10^(-8);
C21new= (C21old)/(kf2*km2);
r21=R21new;
r22=r21;

%MONADA III
R31old=1;
C31old=1;
R32old=Q_3;
C32old=1/Q_3;

kf3=w02;
km3=10^8*C32old/kf3;
R31new=km3*R31old;
R32new=km3*R32old;
C32new=10^(-8);
C31new=(C31old)/(kf3*km3);
r31=R31new;
r32=r31;

%KERDH
k1=1;
k2=2;
k3=2;
H = k1*k2*k3;
aK = 10^(0) / H;
Ra = 100;
Rb = Ra*aK;

%SYNARTHSEIS
T1 = tf( k1 * w00 , [ 0 1 w00 ] );
T2 = tf( k2 * w01^2 , [ 1 ( w01 / Q_2 ) w01^2 ] );
T3 = tf( k3 * w02^2 , [ 1 ( w02 / Q_3 ) w02^2 ] );

TLP = T1*T2*T3;
ltiview ({'bodemag'},TLP)
%before gain regulation
inverseTLP = inv(TLP);
ltiview ({'bodemag'}, inverseTLP)
%after gain regulation
TLP = aK * TLP;
%Attenuation TF
inverseTLP = inv(TLP);
ltiview({'bodemag'}, T1)
ltiview({'bodemag'}, T2)
ltiview({'bodemag'}, T3)
ltiview ({'bodemag'}, TLP)
ltiview({'bodemag'}, T1, T2, T3, TLP)
ltiview ({'bodemag'}, inverseTLP)
plot_transfer_function(TLP , [fp fs] );


%fourier 
Fs= 200*10^3;
T=0.002;
dt=1/Fs;
t=0:dt:(T);
x = sawtooth(2*pi*2000*t);

figure() % άνοιγμα παραθύρου για γραφική παράσταση
plot(t,0.5*(x+1)) % γραφική παράσταση του σήματος

N=T/dt;
xt=lsim(TLP,x,t);
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
y=fft(x,nfft);
figure()
plot(xt)
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*Fs/nfft; 
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

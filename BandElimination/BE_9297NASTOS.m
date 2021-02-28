%AEM: 9297 ON/EP : NASTOS VIKTOR
%BAND ELIMINATION INVERSE CHEBYSHEV
a1=9;
a2=2;
a3=9;
a4=7;

f0=1.8*1000;
f1=1200 + 25*(9-a4);
D= (1/1.8)*(f0^2-f1^2)/f1;
f2= f0^2/f1;
f3= (-D+ sqrt(D^2+4*f0^2))/2;
f4=f0^2/f3;

amin=30-a3;
amax=0.5 + a4/18;

w1 = 2 * pi * f1;
w2= 2 * pi * f2;
w3= 2* pi * f3;
w4 = 2 * pi * f4;
w0 = 2*pi*f0;
Wp = 1;
Ws = ( w2 - w1 ) / ( w4 - w3 );
BW = 2 * pi * ( f2 - f1 );
qc = w0 / BW;

n = acosh( sqrt( ( 10^( amin / 10 ) -1 ) / ( 10^( amax / 10 ) - 1 ) ) ) / acosh( Ws )
n = ceil(n);

e = 1 / sqrt( 10^( amin / 10 ) -1 );
a = ( 1 / n ) * ( asinh( 1 / e ) );
whp = 1 / cosh( acosh( 1 / e ) / n);

ps_k1 = 22.5;
ps_k2 = -22.5;
ps_k3 = 67.5;
ps_k4 = -67.5;

p1 = -sinh(a)*cosd(ps_k1) + j*cosh(a)*sind(ps_k1);
p2 = -sinh(a)*cosd(ps_k2) + j*cosh(a)*sind(ps_k2);
p3 = -sinh(a)*cosd(ps_k3) + j*cosh(a)*sind(ps_k3);
p4 = -sinh(a)*cosd(ps_k4) + j*cosh(a)*sind(ps_k4);

W12 = abs(p1);
W34 = abs(p3);
Q12 = W12/(-2*real(p1));
Q34= W34 /(-2*real(p3));

%antistrofi polwn
InvW12 = 1 / W12;
InvW34 = 1 / W34;
%klimakopoiisi polwn
W12klim=InvW12 * Ws;
W34klim=InvW34 * Ws;
%antistrofi klimakopoiisis polwn
W12klimInv=1/W12klim;
W34klimInv=1/W34klim;
%zeros
wZ1initial = sec((1*pi)/8);
wZ2initial = sec((3*pi)/8);
%klimakopoihsh mhdenikwn
wZ1klim = wZ1initial * Ws;
wZ2klim = wZ2initial * Ws;
%antistrofi polwn
S12= -W12klimInv / (2* Q12);
W12final = sqrt( W12klimInv^2 - S12^2 );
S34= -W34klimInv / (2* Q34);
W34final = sqrt( W34klimInv^2 - S34^2 );
%antistrofi midenikwn
wZ1final= 1 / wZ1klim;
wZ2final= 1 / wZ2klim;

%metasximatismos polou 12
p12Inv = S12 + ( W12 * 1i );
C1= S12^2 + W12final^2;
D1= -2* S12 / qc;
E1= 4 + C1/ qc^2;
G1= sqrt (E1^2 - 4* D1^2);
QS_1= 1/D1 * sqrt ( 1/2* ( E1+ G1) );
k1= -S12 * QS_1 /qc;
W1= k1 + sqrt( k1^2 -1);
w01 = 1/W1 * w0;
w02 = W1 * w0;

%metasximatismos polou 34
p34Inv = S34 + ( W34final * 1i );
C2= S34^2 + W34final^2;
D2= -2* S34 / qc;
E2= 4 + C2/ qc^2;
G2= sqrt (E2^2 - 4* D2^2);
QS_3= 1/D2 * sqrt ( 1/2* ( E2+ G2) );
k2= -S34 * QS_3 /qc;
W3= k2 + sqrt( k2^2 -1);
w03 = 1/W3 * w0;
w04 = W3 * w0;

%metasximatismos fantastikwn midenikwn
Kzero1 = 2 + (wZ1final^2) / (qc^2);
x1 = ( Kzero1 + sqrt( Kzero1^2 - 4 ) ) / 2;
wz11 = w0 * ( sqrt(x1) );
wz21 = w0 / ( sqrt(x1) );

Kzero2 = 2 + (wZ2final^2) / (qc^2);
x2 = ( Kzero2 + sqrt( Kzero2^2 - 4 ) ) / 2;
wz12 = w0 * ( sqrt(x2) );
wz22 = w0 / ( sqrt(x2) );

%Monada I ---> LPN
wz01=wz11/w01;
R11old=1;
C1old=1/(2*QS_1);
R12old=4*QS_1^2;
R15old=(4*QS_1^2)/(wz01^2-1);
koliko1=1/(1 + wz01^2/(2*QS_1^2));
R13old=wz01^2/(2*QS_1^2);
R14old=1;
k1H=1/(R13old+1);
k1L = k1H*(wz01^2);
%klimakopoiisi monadas I
C1new=10^(-7);
km1=C1old/(w01*C1new);
R11new=km1*R11old;
R13new=km1*R13old;
R12new=km1*R12old;
R14new=km1*R14old;
R15new=km1*R15old;

%Monada II ---> HPN
wz02=wz21/w02;
k12=(w02^2/wz21^2) - 1;
k22=((2 + k12)*QS_1^2)/((2+k12)*QS_1^2 + 1);
R21old=1;
R22old=(QS_1^2)*(k12+2)^2;
R23old=1;
R24old=(QS_1^2)*(k12+2);
C22old=1/(QS_1*(2+k12));
C21old=k12*C22old;
k2H=k22*(1/wz02^2);
k2L = k2H*(wz02^2);
%klimakopoiisi monadas II
C22new=10^(-7);
km2=C22old/(w02*C22new);
R21new=km2*R21old;
R23new=km2*R23old;
R22new=km2*R22old;
R24new=km2*R24old;
C21new=C21old/(w02*km2);

%Monada III ---> LPN
wz03=wz12/w03;
R31old=1;
C3old=1/(2*QS_3);
R32old=4*QS_3^2;
R35old=(4*QS_3^2)/(wz03^2-1);
koliko3=1/(1 + wz03^2/(2*QS_3^2));
R33old=wz03^2/(2*QS_3^2);
R34old=1;
k3H=1/(R33old+1);
k3L = k3H*(wz03^2);
%klimakopoiisi monadas III
C3new=10^(-7);
km3=C3old/(w03*C3new);
R31new=km3*R31old;
R33new=km3*R33old;
R32new=km3*R32old;
R34new=km3*R34old;
R35new=km3*R35old;

%Monada IV ---> HPN
wz04=wz22/w04;
k14=(w04^2/wz22^2) - 1;
k24=((2 + k14)*QS_3^2)/((2+k14)*QS_3^2 + 1);
R41old=1;
R42old=(QS_3^2)*(k14+2)^2;
R43old=1;
R44old=(QS_3^2)*(k14+2);
C42old=1/(QS_3*(2+k14));
C41old=k14*C42old;
k4H=k24*(1/wz04^2);
k4L = k4H*(wz04^2);
%klimakopoiisi monadas IV
C42new=10^(-7);
km4=C42old/(w04*C42new);
R41new=km4*R41old;
R43new=km4*R43old;
R42new=km4*R42old;
R44new=km4*R44old;
C41new=C41old/(w04*km4);

%sinartiseis metaforas
T1 = k1H*tf( [1 0 wz11^2], [ 1 (w01/QS_1) w01^2 ] );
T2 = k2H*tf( [1 0 wz21^2], [ 1 (w02/QS_1) w02^2 ] );
T3 = k3H*tf( [1 0 wz12^2], [ 1 (w03/QS_3) w03^2 ] );
T4 = k4H*tf( [1 0 wz22^2], [ 1 (w04/QS_3) w04^2 ] );


%Overall TF
TBE=T1*T2*T3*T4;
%prin ti rithmisi tou kerdous
ltiview ({'bodemag'}, TBE)
%rithmisi kerdous
ktot = k1H*k2H*k3H*k4H;
again = (10^0) / ktot;
TBE = again * TBE
koliko = k1L*k2L*k3L*k4L;
ag=1/koliko;
%ag>1 
Ra = 100;
Rb = Ra*ag;


InverseTBE = inv(TBE);
ltiview({'bodemag'}, T1)
ltiview({'bodemag'}, T2)
ltiview({'bodemag'}, T3)
ltiview({'bodemag'}, T4)
ltiview ({'bodemag'}, TBE)
ltiview({'bodemag'}, T1,T2,T3,T4, TBE)
ltiview ({'bodemag'}, InverseTBE)
plot_transfer_function(TBE, [f1 f2 f3 f4])
plot_transfer_function(InverseTBE, [f1 f2 f3 f4])

f11= (w0- (w0-w3)/3)/(2*pi);
f12=(w0+(w0+w3)/4)/(2*pi);
f13= 0.5*w1 /(2*pi);
f14=(2.4*w2) / (2*pi);
f15= (3*w2)/(2*pi) ;
Fs= 200*10^3;
T=0.002;
dt=1/Fs;
t=0:dt:(T);

u1= cos(2*pi*f11*t)+0.6* cos(2*pi* f12*t)+0.7*cos(2*pi* f13*t)+0.8*cos(2*pi*f14*t)+0.6*cos(2*pi*f15*t);
figure
plot(u1)

N=T/dt;
figure(1)
lsim(TBE,u1,t)
xt=lsim(TBE,u1,t);
figure(2)
plot(t,xt)
n=2^nextpow2(N);
xfourier= fft(xt,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
figure(3)
plot(f,p1)

nfft=n;
y=fft(u1,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*Fs/nfft; 
figure(4)
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




%AEM: 9297 ON/EP : NASTOS VIKTOR
%BAND PASS CHEBYSHEV
a_1 = 9;
a_2 = 2;
a_3 = 9;
a_4 = 7;

f_1 = 400+25*a_3;
f_0 = 650;
f_2 = f_0^2/f_1;
D = 2.3*(f_0^2 -f_1^2)/f_1;
f_3 = (-D+sqrt(D^2+4*f_0^2))/2;
f_4 = f_0^2/f_3;

w_0 = 2*pi*f_0;
w_1 = 2*pi*f_1;
w_2 = 2*pi*f_2;
w_3 = 2*pi*f_3;
w_4 = 2*pi*f_4;

a_min = 27.5+a_4;
a_max = 0.5+(a_3-5)/10;

W_p = 1;
W_s = (w_4 - w_3)/(w_2-w_1);

bw = w_2-w_1;

n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/acosh(W_s);
n = ceil(n);

e = sqrt (10^(a_max/10)-1);

W_hp = cosh(acosh(1/e)/n);

a = asinh(1/e)/n;

ps_k1 = 22.5;
ps_k2 = -22.5;
ps_k3 = 67.5;
ps_k4 = -67.5;

s_k1 = -sinh(a)*cosd(ps_k1) + j*cosh(a)*sind(ps_k1);
s_k2 = -sinh(a)*cosd(ps_k2) + j*cosh(a)*sind(ps_k2);
s_k3 = -sinh(a)*cosd(ps_k3) + j*cosh(a)*sind(ps_k3);
s_k4 = -sinh(a)*cosd(ps_k4) + j*cosh(a)*sind(ps_k4);

Q_1 = sqrt(real(s_k1)^2+imag(s_k1)^2)/abs(2*real(s_k1));
Q_3 = sqrt(real(s_k3)^2+imag(s_k3)^2)/abs(2*real(s_k3));

S_1 = abs(real(s_k1));

qc = w_0/bw;

% Geffe's algorithm
S_1 = abs(real(s_k1));
W_1 = imag(s_k1);
C1 = (S_1)^2+(W_1)^2;
D1 = (2*S_1)/qc;
E1 = 4+C1/(qc)^2;
G1 = sqrt(E1^2 -4*D1^2);
Q_S1 = sqrt((E1+G1)/2)/D1;
k1 = S_1*Q_S1/qc;
W_11 = k1+sqrt(k1^2-1);
w_021 = W_11*w_0;
w_011 =w_0/W_11;

%Poloi meta ton metasximatismo
sk1 = w_011 / (2*Q_S1);
sk2 = w_021 / (2*Q_S1);
wk1 = w_011 * sin(acos(1/(2*Q_S1)));
wk2 = w_021 * sin(acos(1/(2*Q_S1)));
yk1 = acos(1/(2*Q_S1))*180/pi;
yk2 = -acos(1/(2*Q_S1))*180/pi;

% Geffe's algorithm
S_3 = abs(real(s_k3));
W_3 = imag(s_k3);
C3 = (S_3)^2+(W_3)^2;
D3 = 2*S_3/qc;
E3 = 4+C3/(qc)^2;
G3 = sqrt(E3^2 -4*D3^2);
Q_S3 = sqrt((E3+G3)/2)/D3;
k3 = S_3*Q_S3/qc;
W_33 = k3+sqrt(k3^2-1);
w_023 = W_33*w_0;
w_013 =w_0/W_33;

%Poloi meta ton metasximatismo
sk3 = w_013 / (2*Q_S3);
sk4 = w_023 / (2*Q_S3);
wk3 = w_013 * sin(acos(1/(2*Q_S3)));
wk4 = w_023 * sin(acos(1/(2*Q_S3)));
yk3 = acos(1/(2*Q_S3))*180/pi;
yk4 = -acos(1/(2*Q_S3))*180/pi;

%MONADA I (Q enhancement) Q>>5
C_11 = 1;
C_12 = C_11;
b1 = 1;
R_11 = 1/sqrt(b1);
R_12 = sqrt(b1);
k_1 = (Q_S1*(b1+2)-sqrt(b1))/(2*Q_S1-sqrt(b1));
R_A1 = 1;
R_B1 = (k_1-1)*R_A1;
%KLIMAKOPOIHSH MONADAS 1
kf_1 = w_011;
C_n11 = 10^(-7);
km_1 = C_11/(C_n11*kf_1);
C_n12 = C_n11;
R_n11 = km_1 * R_11;
R_n12 = km_1 * R_12;
R_A1=km_1*R_A1;
R_B1=km_1*R_B1;


%MONADA II
C_21 = 1;
C_22 = C_21;
b2 = 1;
R_21 = 1/sqrt(b2);
R_22 = sqrt(b2);
k_2 = (Q_S1*(b2+2)-sqrt(b2))/(2*Q_S1-sqrt(b2));
R_A2 = 1;
R_B2 = (k_2-1)*R_A2;
%KLIMAKOPOIHSH
kf_2 = w_021;
C_n21 = 10^(-7);
km_2 = C_21/(C_n21*kf_2);
C_n22 = C_n21;
R_n21 = km_2 * R_21;
R_n22 = km_2 * R_22;
R_A2=km_2*R_A2;
R_B2=km_2*R_B2;

%MONADA III (Q enhancement) Q>>5
C_31 = 1;
C_32 = C_31;
b3 = 1;
R_31 = 1/sqrt(b3);
R_32 = sqrt(b3);
k_3 = (Q_S3*(b3+2)-sqrt(b3))/(2*Q_S3-sqrt(b3));
R_A3 = 1;
R_B3 = (k_3-1)*R_A3;
%KLIMAKOPOIHSH MONADAS III
kf_3 = w_013;
C_n31 = 10^(-7);
km_3 = C_31/(C_n31*kf_3);
C_n32 = C_n31;
R_n31 = km_3 * R_31;
R_n32 = km_3 * R_32;
R_A3=km_3*R_A3;
R_B3=km_3*R_B3;

%MONADA IV
C_41 = 1;
C_42 = C_41;
b4 = 1;
R_41 = 1/sqrt(b4);
R_42 = sqrt(b4);
k_4 = (Q_S3*(b4+2)-sqrt(b4))/(2*Q_S3-sqrt(b4));
R_A4 = 1;
R_B4 = (k_4-1)*R_A4;
%KLIMAKOPOIHSH MONADAS IV
kf_4 = w_023;
C_n41 = 10^(-7);
km_4 = C_41/(C_n41*kf_4);
C_n42 = C_n41;
R_n41 = km_4 * R_41;
R_n42 = km_4 * R_42;
R_A4=km_4*R_A4;
R_B4=km_4*R_B4;

H_1 = k_1*b1/(2*(k_1-1)-b1);
H_2 = k_2*b2/(2*(k_2-1)-b2);
H_3 = k_3*b3/(2*(k_3-1)-b3);
H_4 = k_4*b4/(2*(k_4-1)-b4);

%SYNARTHSEIS METAFORAS
TF_1 = tf ([0 H_1*(w_011/Q_S1) 0], [1 w_011/Q_S1 w_011^2]);
TF_2 = tf ([0 H_2*(w_021/Q_S1) 0], [1 w_021/Q_S1 w_021^2]);
TF_3 = tf ([0 H_3*(w_013/Q_S3) 0], [1 w_013/Q_S3 w_013^2]);
TF_4 = tf ([0 H_4*(w_023/Q_S3) 0], [1 w_023/Q_S3 w_023^2]);

H_21 = bode(TF_1,w_0);
H_22 = bode(TF_2,w_0);
H_23 = bode(TF_3,w_0);
H_24 = bode(TF_4,w_0);

%RYTHMISH KERDOYS

Holiko=bode(TF_1,w_0)*bode(TF_2,w_0)*bode(TF_3,w_0)*bode(TF_4,w_0);
aK = 10^(0.25)/Holiko;
za = R_n11/aK;
zb = R_n11/(1-aK);

%ltiview ({'bodemag'}, TBP)
%inverseTBP = inv(TBP);
%ltiview ({'bodemag'}, inverseTBP)
%AFTER GAIN REGULATION
 %TBP = aK * TBP;
 
 TBP = TF_1*TF_2*TF_3*TF_4;
 plot_transfer_function(TBP , [f_1 f_2 f_3 f_4] );
 ltiview({'bodemag'}, TF_1)
 ltiview({'bodemag'}, TF_2)
 ltiview({'bodemag'}, TF_3)
 ltiview({'bodemag'}, TF_4)
 ltiview ({'bodemag'}, TBP);
 ltiview({'bodemag'}, TF_1, TF_2, TF_3, TF_4, TBP)
 InvSys_new = inv (TBP);
 ltiview ({'bodemag'}, 1/TBP)
 plot_transfer_function(1/TBP , [f_1 f_2 f_3 f_4] );
 TBP = aK*TF_1*TF_2*TF_3*TF_4;
 plot_transfer_function(1/TBP , [f_1 f_2 f_3 f_4] );
 plot_transfer_function(TBP , [f_1 f_2 f_3 f_4] );
 ltiview ({'bodemag'}, TBP)
 %plot_transfer_function(InvSys_new , [f_1 f_2] );


%eisodos k fft
multisim_w0 = w_0 - ((w_0 - w_1)/3);
multisim_w1 = w_0 + ((w_0 + w_1)/4);
multisim_w2 = 0.5 * w_3;
multisim_w3 = 2.4 * w_4;
multisim_w4 = 3 * w_4;
syms t;
input_signal(t) = cos(multisim_w0*t) + 0.6*cos(multisim_w1*t) + 0.7*cos(multisim_w2*t) ...
               + 0.8*cos(multisim_w3*t) + 0.6*cos(multisim_w4*t);
figure;
time_sym = 0:0.001:0.15;
ezplot(input_signal(t),time_sym);
Fs = 20000;          % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;            % Length of signal
time = (0:L-1)*T;     % Time vector           

x = cos(multisim_w0*time ) + 0.6*cos(multisim_w1*time ) + cos(multisim_w2*time ) ...
               + 0.8*cos(multisim_w3*time ) + 0.4*cos(multisim_w4*time );

Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

sim = lsim(TBP,x,time);
Y = fft(sim);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
multisim_f0 = multisim_w0 / (2*pi);
multisim_f1 = multisim_w1 / (2*pi);
multisim_f2 = multisim_w2 / (2*pi);
multisim_f3 = multisim_w3 / (2*pi);
multisim_f4 = multisim_w4 / (2*pi);

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

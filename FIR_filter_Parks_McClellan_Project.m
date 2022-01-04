%% Design FIR filter using Parks-McClellan algorithm
clear;
clc;
disp('#---------------- H1_Filter_Result ----------------#');
% ========================================================== %
%                       -- H1 Filter --
% ========================================================== %
%% Low-pass filter H1 Specify
fpass_edge_h1 = 20e3;                  % Passband edge frequency (kHz)
fstop_edge_h1 = 40e3;                  % Stopband edge frequency (kHz)
rate = 250e3;                          % Sampling rate (kHz)
Rpass_H1 = 0.4;                           % Passband ripple (dB) 
Rstop_H1 = 80;                            % Stopband ripple (dB)
xp = fpass_edge_h1/rate;
xs = fstop_edge_h1/rate;
w_pass1 = (2*pi) * xp;
w_stop1 = (2*pi) * xs;
% Compute deviations
dev_H1 = [((10^(Rpass_H1/20)-1)/(10^(Rpass_H1/20)+1))/5, 10^(-Rstop_H1/20)/7];  

%% Low-pass filter H1 design
[M, f0, A0, W] = firpmord([fpass_edge_h1 fstop_edge_h1], [1 0], dev_H1, rate);
% firpmord(f, a, dev, fs):      f:[fpass_edge, fstop_edge];   a:[1, 0]:low_pass_filter;
% dev(ripple_error):[(10^(rp/20)-1)/(10^(rp/20)+1), 10^(-rs/20)];     fs:sampling_rate;

disp('# H1 filter order (M):');
disp(M);   % M: order
h1 = firpm(M, f0, A0, W);
% fvtool(h,1)
% firpm(n, f0, a0, w):      n:order;      f0: each Truncated-edge frequency ; 
% a0: each Truncated-edge frequency corresponding amplitudes;   
% w: weights (強調 minimum 某些頻段相對於其他頻段的誤差, Stopband ripple 比 Passband ripple 少 10 倍) -> 壓縮 stopband ripple 讓它變更小
H1 = 20 * log10(abs(fft(h1, 10000)));
H1_gain_dB = 20 * log10(2);
H1 = (H1 - max(H1)) + H1_gain_dB ;
f = [0:4999]/10000;

figure(1);
plot(f, H1(1:5000), [xs 0.4], [-80 -80], 'k:', [xs xs], [-10 -80], 'k:');
axis([0 0.3, -120 10]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('H1 Filter designed by Parks-McClellan"s method');

figure(2);
plot(f, H1(1:5000), [0 xp], [5.6 5.6], 'k:', [0 xp], [6.4 6.4], 'k:', [xp xp], [5.6 -1], 'k:', [xp xp], [0 1], 'k:');
axis([0 xp*1.2, 4.5 6.6]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('Passband');

%% ---- H1 fixed-point approximation ( Quantized )
n = 15;
% h_fix = floor(h1/max(h1) * 2^n + 0.5);
h1_fix = round(h1/max(h1) * 2^n) / 2^n;
H1fix = abs(fft(h1_fix, 10000));
H1_fix = 20 * log10(H1fix);
H1_fix = (H1_fix - max(H1_fix)) + H1_gain_dB ;

figure(3);
plot(f, H1_fix(1:5000), [xs 0.4], [-80 -80], 'k:', [xs xs], [-10 -80], 'k:');
axis([0 0.3 -120 10]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('H1 Filter 15-bits-fixed-point version');

figure(4);
plot(f, H1_fix(1:5000), [0 xp], [5.6 5.6], 'k:', [0 xp], [6.4 6.4], 'k:', [xp xp], [5.6 -1], 'k:', [xp xp], [0 1], 'k:');
axis([0 xp*1.2 4.5 6.6]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('Passband (fixed-point version)');

h1_fix_int = h1_fix * 2^n;    % filter-coefficients = filter-total-length = order + 1 ;  % order = delay-register (Z^-1) number
disp('# H1_fixed_point_Quantized:');
disp(h1_fix_int);


%% H1 add Group Delay
filter_length = M + 1;                            % M (Order) must be odd number ; M (order) +1 = filter_length
Delta = zeros(1,filter_length);                   % Coefficient add a delta vector = filter add Group delay
Delta(1,(filter_length +1)/2) = 1;                % Center point of symmetry coefficient set to 1, else set to 0
% grpdelay(B, A, n, fs):  H = B/A (transfer_function, when A = 1, filter = FIR)
% n = points to be response; 
h1_mix_delay = h1_fix_int + Delta;                % H1_filter_coefficients add H1_group_delays
disp('# H1_fixed_point_Add_Delay:');
disp(h1_mix_delay);

H1_mix_delay = 20 * log10(abs(fft(h1_mix_delay, 10000)));
H1_mix_delay = (H1_mix_delay - max(H1_mix_delay)) + H1_gain_dB ;
f = [0:4999]/10000;

figure(5);
plot(f, H1_mix_delay(1:5000), [xs 0.4], [-80 -80], 'k:', [xs xs], [-10 -80], 'k:');
axis([0 0.3 -120 10]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('H1 Filter designed by Group Delay');

figure(6);
plot(f, H1_mix_delay(1:5000), [0 xp], [5.6 5.6], 'k:', [0 xp], [6.4 6.4], 'k:', [xp xp], [5.6 -1], 'k:', [xp xp], [0 1], 'k:');
axis([0 xp*1.2 4.5 6.6]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('Passband');

%%
% Finish H1 & H1 delay filter
%%
disp('#---------------- H2_Filter_Result ----------------#');
% ========================================================== %
%                       -- H2 Filter --
% ========================================================== %
%% Low-pass filter H2 Specify
fpass_edge_h2 = 60e3;                  % Passband edge frequency (kHz)
fstop_edge_h2 = 80e3;                  % Stopband edge frequency (kHz)
rate = 250e3;                          % Sampling rate (kHz)
Rpass_H2 = 0.4;                        % Passband ripple (dB) 
Rstop_H2 = 80;                         % Stopband ripple (dB)
xp = fpass_edge_h2/rate;
xs = fstop_edge_h2/rate;
w_pass2 = (2*pi) * xp;
w_stop2 = (2*pi) * xs;
% Compute deviations
dev_H2 = [((10^(Rpass_H1/20)-1)/(10^(Rpass_H1/20)+1))/2, 10^(-Rstop_H1/20)/3];

%% Low-pass filter H2 design
[M, f0, A0, W] = firpmord([fpass_edge_h2 fstop_edge_h2], [1 0], dev_H2, rate);
% firpmord(f, a, dev, fs):      f:[fpass_edge, fstop_edge];   a:[1, 0]:low_pass_filter;
% dev(ripple_error):[(10^(rp/20)-1)/(10^(rp/20)+1), 10^(-rs/20)];     fs:sampling_rate;

disp('# H2 filter order (M):');
disp(M);   % M: order
h2 = firpm(M, f0, A0, W);
% fvtool(h,1)
% firpm(n, f0, a0, w):      n:order;      f0: each Truncated-edge frequency ; 
% a0: each Truncated-edge frequency corresponding amplitudes;   
% w: weights (強調 minimum 某些頻段相對於其他頻段的誤差, Stopband ripple 比 Passband ripple 少 10 倍) -> 壓縮 stopband ripple 讓它變更小
H2 = 20 * log10(abs(fft(h2, 10000)));
H2_gain_dB = 20 * log10(1);
H2 = (H2 - max(H2)) + H2_gain_dB ;
f = [0:4999]/10000;

figure(7);
plot(f, H2(1:5000), [xs 0.4], [-80 -80], 'k:', [xs xs], [-10 -80], 'k:');
axis([0 0.3 -120 10]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('H2 Filter designed by Parks-McClellan"s method');

figure(8);
plot(f, H2(1:5000), [0 xp], [-0.4 -0.4], 'k:', [0 xp], [0 0], 'k:', [xp xp], [-0.4 -1], 'k:', [xp xp], [0 1], 'k:');
axis([0 xp*1.2 -1.2 0.2]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('Passband');

%% ---- H2 fixed-point approximation ( Quantized )
n = 15;
% h_fix = floor(h/max(h) * 2^n + 0.5);
h2_fix = round(h2/max(h2) * 2^n) / 2^n;
H2fix = abs(fft(h2_fix, 10000));
H2_fix = 20 * log10(H2fix);
H2_fix = (H2_fix - max(H2_fix)) + H2_gain_dB ;

figure(9);
plot(f, H2_fix(1:5000), [xs 0.4], [-80 -80], 'k:', [xs xs], [-10 -80], 'k:');
axis([0 0.3 -120 10]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('H2 Filter 15-bits-fixed-point version');

figure(10);
plot(f, H2_fix(1:5000), [0 xp], [-0.4 -0.4], 'k:', [0 xp], [0 0], 'k:', [xp xp], [-0.4 -1], 'k:', [xp xp], [0 1], 'k:');
axis([0 xp*1.2 -1.2 0.2]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('Passband (fixed-point version)');

h2_fix_int = h2_fix * 2^n;    % filter-coefficients = filter-total-length = order + 1 ;  % order = delay-register (Z^-1) number
disp('# H2_fixed_point_Quantized:');
disp(h2_fix_int);

%% Mix two filter
system_filter = H1_fix + H2_fix ;

figure(11);
plot(f, system_filter(1:5000), [xs 0.4], [-160 -160], 'k:', [xs xs], [-10 -160], 'k:', [xs 0.16], [-80 -80], 'k:', [xs xs], [0 -80], 'k:');
axis([0 0.5 -300 10]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('H1 + H2 Filter 15-bits-fixed-point version');

figure(12);
plot(f, system_filter(1:5000), [0 xp], [5.6 5.6], 'k:', [0 xp], [6.4 6.4], 'k:', [xp xp], [5.6 -1], 'k:', [xp xp], [0 1], 'k:');
axis([0 xp*1.2 4.5 6.6]);
xlabel('normalized frequency (cycle/sample)');
ylabel('normalized magnitude response (dB)');
title('Passband (fixed-point version)');





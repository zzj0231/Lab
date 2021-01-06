function Hd = filter8M_start1500k_stop1700k
%FILTER8M_START1500K_STOP1700K 返回离散时间滤波器对象。

% MATLAB Code
% Generated by MATLAB(R) 9.8 and DSP System Toolbox 9.10.
% Generated on: 19-Dec-2020 11:17:19

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in MHz.
Fs = 8.125;  % Sampling Frequency

Fpass = 1.5;             % Passband Frequency
Fstop = 1.7;             % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.001;           % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
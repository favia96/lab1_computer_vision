%% LAB 1 - COMPUTER VISION
%% by Federico Favia, Martin De Pellegrini

clear
clc
clear all

%% Part 1: Properties of discrete Fourier Transform

% Fhat = zeros(128, 128);
% p = 3;
% q = 2;
% Fhat(p , q ) = 1;
% 
% F = ifft2(Fhat);
% 
% Fabsmax = max(abs(F(:)));
% showgrey(real(F), 64, -Fabsmax, Fabsmax)
% showgrey(imag(F), 64, -Fabsmax, Fabsmax)
% showgrey(abs(F), 64, -Fabsmax, Fabsmax)
% showgrey(angle(F), 64, -pi, pi)

[Fhat, F, Fabsmax] = fftwave(128,128);

%% Part 2: Gaussian convolution implemented via FFT

%% Part 3: Smoothing

% office = office256;
% add = gaussnoise(office, 16);
% sap = sapnoise(office, 0.1, 255);
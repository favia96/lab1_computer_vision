%% LAB 1 - COMPUTER VISION, November 2019
%% by Federico Favia, Martin De Pellegrini

%% Initialization
clear ; close all; clc

%% Part 1: Properties of discrete Fourier Transform
% basis functions

p = [1, 2, 5, 9, 17, 17, 5, 64, 120, 70, 125, 128];
q = [1, 2, 9, 5, 9, 121, 1, 64, 70, 120, 1, 128];

for i = 1 : length(p)
    figure()
    [Fhat_, F_, Fabsmax_, amplitude_, wavelength_ ] = fftwave(p(i),q(i));
    w = waitforbuttonpress
end
    

%% linearity

% F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)]; %test image
% G = F'; %transposed
% H = F + 2 * G;
% 
% for i = 1 : 3
%     figure()
%     switch i
%         case 1
%             showgrey(F);
%             Fhat = fft2(F);
%             figure()
%             showgrey(1 + abs(Fhat));
%         case 2
%             showgrey(G);
%             Ghat = fft2(G);
%             figure()
%             showgrey(log(1 + abs(Ghat)));
%         otherwise
%             showgrey(H);
%             Hhat = fft2(H);
%             figure()
%             %showgrey(log(1 + abs(Hhat)));
%             showgrey(log(1 + abs(fftshift(Hhat))));
%     end
% end

%% multiplication

%% Part 2: Gaussian convolution implemented via FFT

%% Part 3: Smoothing

% office = office256;
% add = gaussnoise(office, 16);
% sap = sapnoise(office, 0.1, 255);
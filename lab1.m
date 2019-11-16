%% LAB 1 - COMPUTER VISION, November 2019
%% by Federico Favia, Martin De Pellegrini

%% Initialization
clear ; close all; clc

%% Part 1: Properties of discrete Fourier Transform
% basis functions

% p = [1, 2, 5, 9, 17, 17, 5, 64, 120, 70, 125, 128];
% q = [1, 2, 9, 5, 9, 121, 1, 64, 70, 120, 1, 128];
% 
% for i = 1 : length(p)
%     figure()
%     [Fhat_, F_, Fabsmax_, amplitude_, wavelength_ ] = fftwave(p(i),q(i));
%     w = waitforbuttonpress
% end
%     

%% linearity
%%F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)]; %test image
%%G = F'; %transposed
% H = F + 2 * G;
% 
% figure();
% subplot(1,4,1);
% imshow(F);
% title("Original F");
%%Fhat = fft2(F);
% subplot(1,4,2);
% showgrey(abs(Fhat));
% title("Spectra F");
% subplot(1,4,3);
% showgrey(log(1 + abs(Fhat)));
% title("Log spectra F");
% subplot(1,4,4);
% showgrey(log(1 + abs(fftshift(Fhat))));
% title("Log spectra shift F");
% 
% figure();
% subplot(1,4,1);
% showgrey(G);
% title("Transposed F = G");
%%Ghat = fft2(G);
% subplot(1,4,2);
% showgrey(abs(Ghat));
% title("Spectra G");
% subplot(1,4,3);
% showgrey(log(1 + abs(Ghat)));
% title("Log spectra G");
% subplot(1,4,4);
% showgrey(log(1 + abs(fftshift(Ghat))));
% title("Log spectra shift G");
% 
% figure();
% subplot(1,5,1);
% showgrey(H);
% title("Original H");
% Hhat = fft2(H);
% subplot(1,5,2);
% showgrey(abs(Hhat));
% title("Spectra G");
% subplot(1,5,3);
% showgrey(log(1 + abs(Hhat)));
% title("Log spectra H");
% subplot(1,5,4);
% showgrey(log(1 + abs(fftshift(Hhat))));
% title("Log spectra shift H");
% subplot(1,5,5);
% showgrey(ifft2(Fhat+2*Ghat));
% title("linearity");

%% multiplication
% figure();
% subplot(1,5,1);
% Z = F .* G;
% Zhat = fft2(Z);
% showgrey(Z);
% title("f x g");
% subplot(1,5,2);
% showgrey(abs(Zhat));
% title("spectra f x g");
% subplot(1,5,3);
% showgrey(log(1+abs(Zhat)));
% title("log spectra f x g");
% subplot(1,5,4);
% showgrey(log(1+abs(fftshift(Zhat))));
% title("log spectra shift f x g");
% subplot(1,5,5);
% showfs(Zhat);
% title("fs f x g");
% 
% Zhat2=conv2(Fhat,Ghat,'same');
% figure();
% Z2 = ifft2(Zhat);
% subplot(1,3,1)
% showgrey(Z2); title("z2=f x g");
% subplot(1,3,2)
% title("Zhat2");
% subplot(1,3,3)
% Zhat3=fft2(Z2);
% showfs(Zhat3); title("fft2(Z2)");

%% information in fourier phase and magnitude
img = phonecalc128;
figure;
showgrey(img);
a=(10.0)^-10;
figure;
showgrey(pow2image(img,a));
randphaseimage(;

%% Part 2: Gaussian convolution implemented via FFT

%% Part 3: Smoothing

% office = office256;
% add = gaussnoise(office, 16);
% sap = sapnoise(office, 0.1, 255);
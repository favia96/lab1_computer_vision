%% LAB 1 - COMPUTER VISION
%% by Federico Favia, Martin De Pellegrini

clear
clc
clear all

%% Part 1: Properties of discrete Fourier Transform

p = [1, 5, 9, 17, 17, 5, 64, 125, 128];
q = [1, 9, 5, 9, 121, 1, 64, 1, 128];


% for i = 1 : length(p)
%   	Fhat = zeros(128, 128);
%     Fhat(p(i) , q(i) ) = 1;
%     F = ifft2(Fhat);
%     %figure()
%     %set(FigH,sprintf('(p, q) = (%d, %d)', p(i), q(i)))
%     Fabsmax = max(abs(F(:)));
%     
%     subplot(2,2,1)
%     showgrey(real(F), 64, -Fabsmax, Fabsmax)
%         
%     subplot(2,2,2)
%     showgrey(imag(F), 64, -Fabsmax, Fabsmax)
%     
%     subplot(2,2,3)
%     showgrey(abs(F), 64, -Fabsmax, Fabsmax)
%     
%     subplot(2,2,4)
%     showgrey(angle(F), 64, -pi, pi)
%     
%     keyboard
% end


for i = 1 : length(p)
    figure()
    [real_, Fhat_, F_, Fabsmax_, amplitude_, wavelength_ ] = fftwave(p(i),q(i));
    keyboard
end
    


%% Part 2: Gaussian convolution implemented via FFT

%% Part 3: Smoothing

% office = office256;
% add = gaussnoise(office, 16);
% sap = sapnoise(office, 0.1, 255);
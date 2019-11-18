%% LAB 1 - COMPUTER VISION
%% by Federico Favia, Martin De Pellegrini

clear
clc
% clear all

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

% p = [5 9 17 17 5 125];
% q = [9 5 9 121 1 1];

p = [1 64 25]; q = [1 64 95];
for i = 1:length(p)
    figure;
    [Fhat, F, Fabsmax] = fftwave(p(i),q(i));
end
% [Fhat, F, Fabsmax] = fftwave(50,43);

%% question 10
F = [zeros(56, 128); ones(16, 128); zeros(56,128)];  %showgray(F);
G = F';

figure;
subplot(121);
showgray(F.*G);
subplot(122);
showfs(fft2(F.*G));
X = fft2(F.*G);
figure; showfs(X);

%the result we obtain from the moltiplication is a quare in the center of
%the image, the relative fourier transform is a sinc in the frequency
%domain.

%it is possible to achieve the same result by first computing the fft of
%the two images separately, then compute the convolution between the two.

%Q10
Fhat = fft2(F); Ghat = fft2(G);
figure;
subplot(121); showfs(Fhat);
subplot(122); showfs(Ghat);

H = conv2(Fhat,Ghat,'same');
figure;
showfs(H);
% figure;
% showgray(log(1+abs(H)));

%% Q11
F1 = [zeros(60,128); ones(8,128); zeros(60,128)].*[zeros(128,48) ones(128,32) zeros(128,48)];
figure; showgray(F1);
figure; showfs(fft2(F1));
figure;
subplot(211); showfs(fft2(F)); title('question 10');
subplot(212); showfs(fft2(F1)); title('question 11');

%% Q12 rotation
alpha = [30 45 60 90];
% figure;
% subplot(221); showgray(F1);
% subplot(222); showfs(fft2(F1));

G1 = rot(F1, alpha(1));
Ghat1 = fft2(G1);
G2 = rot(F1, alpha(2));
Ghat2 = fft2(G2);
G3 = rot(F1, alpha(3));
Ghat3 = fft2(G3);
G4 = rot(F1, alpha(4));
Ghat4 = fft2(G4);


subplot(211); showgray(G1); axis on
subplot(212); showfs(Ghat);

Hhat1 = rot(fftshift(Ghat1), -alpha(1));
Hhat2 = rot(fftshift(Ghat2), -alpha(2));
Hhat3 = rot(fftshift(Ghat3), -alpha(3));
Hhat4 = rot(fftshift(Ghat4), -alpha(4));

figure; 
subplot(241); showfs(fft2(F1)); title('original not rotate');
subplot(242); showfs(fft2(F1)); title('original not rotate');
subplot(243); showfs(fft2(F1)); title('original not rotate');
subplot(244); showfs(fft2(F1)); title('original not rotate');
subplot(245); showgray(log(1+abs(Hhat1))); title('rotated back');
subplot(246); showgray(log(1+abs(Hhat2))); title('rotated back');
subplot(247); showgray(log(1+abs(Hhat3))); title('rotated back');
subplot(248); showgray(log(1+abs(Hhat4))); title('rotated back');

%% Q13
im = phonecalc128;
a = 10.0^(-10);
figure;
subplot(121); showgray(pow2image(im,a));
subplot(122); showgray(randphaseimage(im));

%% Part 2: Gaussian convolution implemented via FFT
%write on e matlab procedure that perform a gaussian filtering in fourier
%domain
var = [0.1 0.3 1.0 10.0 100.0];

im = phonecalc128;
% im_out = gaussfft(im, 1.0);
[im_out1, varMat1] = gaussfft(im, var(1));
[im_out2, varMat2] = gaussfft(im, var(2));
[im_out3, varMat3] = gaussfft(im, var(3));
[im_out4, varMat4] = gaussfft(im, var(4));
[im_out5, varMat5] = gaussfft(im, var(5));

figure;
subplot(321); showgray(im); title('Original image');
subplot(322); showgray(im_out1); title('Variance in frequency = 0.1');
subplot(323); showgray(im_out2); title('Variance in frequency = 0.3');
subplot(324); showgray(im_out3); title('Variance in frequency = 1.0');
subplot(325); showgray(im_out4); title('Variance in frequency = 10.0');
subplot(326); showgray(im_out5); title('Variance in frequency = 100.0');



%% test
im = phonecalc128;
var = 100.0;
Im_in = fft2(im);  %compute the fourier transform
    
[dim1,dim2] = size(Im_in);
    
% [x,y] = meshgrid(0:dim1-1, 0:dim2-1);  %mesh for the Gaussian Kernel
x = -dim1/2:dim1/2 -1;
y = -dim2/2:dim2/2 -1;
    
 %Gaussian Kernel
gauss = (1/2*pi*var).*exp((-(x.^2  + y.^2)./(2*var)));
% gauss = fftshift(gauss);

% gauss = (1/2*pi*var).*exp(-(x.^2)/(2*var)).*exp(-(y.^2)/(2*var));
% gauss = (1/2*pi*var).*exp((-((x-dim1/2).^2  + (y-dim2/2).^2)./(2*var)));
% gauss = fftshift(gauss);

mesh(gauss);


%filtering
Im_out = Im_in.*gauss;
    
%output
im_out = ifft2(Im_out);

figure;
subplot(121); showgray(im);
subplot(122); showgray(im_out);

%% Part 3: Smoothing

% office = office256;
% add = gaussnoise(office, 16);
% sap = sapnoise(office, 0.1, 255);
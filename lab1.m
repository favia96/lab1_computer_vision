%% LAB 1 - COMPUTER VISION, November 2019
%% by Federico Favia, Martin De Pellegrini

%% Initialization
clear ; close all; clc

%% Part 1: Properties of discrete Fourier Transform (Q1-Q2-Q3-Q4-Q5)
% basis functions

% p = [1, 2, 10, 45, 80, 17, 17, 5, 64, 10, 70, 125, 128];
% q = [1, 2, 80, 45, 10, 9, 121, 1, 64, 70, 120, 1, 128];
% 
% for i = 1 : length(p)
%      figure()
%      [Fhat_, F_, Fabsmax_, amplitude_, wavelength_ ] = fftwave(p(i),q(i));
%      w = waitforbuttonpress
% end

% Q2
Q = 128; % dimension of original image
orig_image = zeros(Q,Q);
period = 5; %number of periods
for m = 1 : 1 : Q % sine wave in spatial domain (original image)
    for n = 1 : 1 : Q
        orig_image(m,n) = 0.5 + 0.5.*sin(2.*pi.*(n./(Q./period))); %the "0.5" constants give values between 0 and 1 for use with "imshow"
    end
end
subplot(1,2,1); imshow(orig_image); title('Original image');
B = abs(fftshift(fft2(orig_image-0.5))); %delete the added scaling values before FFT, gives values symmetric around the x_axis
subplot(1,2,2); imshow(B); title('FFT2 Magnitude');  

%% linearity (Q7-Q8-Q9)
F = [ zeros(56, 128); ones(16, 128); zeros(56, 128)]; %test image
G = F'; %transposed
% H = F + 2 * G;
% 
% figure();
% subplot(1,4,1); showgrey(F); title("Original F");
Fhat = fft2(F);
% subplot(1,4,2); showgrey(abs(Fhat)); title("Spectra F");
% subplot(1,4,3); showgrey(log(1 + abs(Fhat))); title("Log spectra F");
% subplot(1,4,4); showgrey(log(1 + abs(fftshift(Fhat)))); title("Log spectra shift F");
% figure();
% subplot(1,4,1); showgrey(G); title("Transposed F = G");
Ghat = fft2(G);
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
% subplot(1,5,1); showgrey(H); title("Original H");
% Hhat = fft2(H);
% subplot(1,5,2); showgrey(abs(Hhat)); title("Spectra G");
% subplot(1,5,3); showgrey(log(1 + abs(Hhat))); title("Log spectra H");
% subplot(1,5,4); showgrey(log(1 + abs(fftshift(Hhat)))); title("Log spectra shift H");
% subplot(1,5,5); showgrey(ifft2(Fhat+2*Ghat)); title("linearity");

%% multiplication (Q10)
% Z = F .* G;
% Zhat = fft2(Z); 
% figure();
% subplot(1,5,1); showgrey(Z); title("f x g");
% subplot(1,5,2); showgrey(abs(Zhat)); title("spectra f x g");
% subplot(1,5,3); showgrey(log(1+abs(Zhat))); title("log spectra f x g");
% subplot(1,5,4); showgrey(log(1+abs(fftshift(Zhat)))); title("log spectra shift f x g");
% subplot(1,5,5); showfs(Zhat); title("fs f x g");
% 
% Zhat2 = conv2(fftshift(Fhat),fftshift(Ghat),'same');
% Z2 = ifft2(Zhat); % approximation of Z = F. * G
% figure();
% subplot(1,3,1); % showgrey(Z2); title("z2=f x g");
% test2 = abs(fftshift(Zhat2));
% subplot(1,3,2); showgrey(log(1+abs(fftshift(test2)))); %* problems title("Zhat2");
% Zhat3 = fft2(Z2);
% subplot(1,3,3); showfs(Zhat3); title("fft2(Z2)");

%% scaling (Q11)
F1 = [zeros(60,128); ones(8,128); zeros(60,128)] .* [zeros(128,48) ones(128,32) zeros(128,48)];
F1hat = fft2(F1);
% figure()
% subplot(2,2,1); showgrey(F1); title('image F1');
% subplot(2,2,2); showfs(F1hat); title('log spectrum shifted F1')
% subplot(2,2,3); showgrey(F); title('image F');
% subplot(2,2,4); showfs(Fhat); title('log spectrum shifted F');

%% rotation (Q12)
% alpha = [30, 45, 60, 90];
% % figure;
% % subplot(221); showgrey(F1);
% % subplot(222); showfs(fft2(F1));

% G1 = rot(F1, alpha(1));
% Ghat1 = fft2(G1);
% G2 = rot(F1, alpha(2));
% Ghat2 = fft2(G2);
% G3 = rot(F1, alpha(3));
% Ghat3 = fft2(G3);
% G4 = rot(F1, alpha(4));
% Ghat4 = fft2(G4);
% 
% subplot(211); showgrey(G1); axis on
% subplot(212); showfs(Ghat);
% 
% Hhat1 = rot(fftshift(Ghat1), -alpha(1));
% Hhat2 = rot(fftshift(Ghat2), -alpha(2));
% Hhat3 = rot(fftshift(Ghat3), -alpha(3));
% Hhat4 = rot(fftshift(Ghat4), -alpha(4));
% 
% figure; 
% subplot(241); showfs(fft2(F1)); title('original not rotate');
% subplot(242); showfs(fft2(F1)); title('original not rotate');
% subplot(243); showfs(fft2(F1)); title('original not rotate');
% subplot(244); showfs(fft2(F1)); title('original not rotate');
% subplot(245); showgrey(log(1+abs(Hhat1))); title('rotated back');
% subplot(246); showgrey(log(1+abs(Hhat2))); title('rotated back');
% subplot(247); showgrey(log(1+abs(Hhat3))); title('rotated back');
% subplot(248); showgrey(log(1+abs(Hhat4))); title('rotated back');

% %% Information in fourier phase and magnitude (Q13)
% img = nallo128; % try with phonecalc128, few128, nallo128
% a = 10.0^(-10);
% figure()
% subplot(1,3,1); showgrey(img); title('Original image');
% subplot(1,3,2); showgrey(pow2image(img,a)); title('Phase'); % phase preserved and magnitude replaced by non-lin fct
% subplot(1,3,3); showgrey(randphaseimage(img)); title('Magnitude'); % magnitude kept and phase as random distr

%% Part 2: Gaussian convolution implemented via FFT (Q14-Q15-Q16)
vars = [0.1 0.3 1.0 10.0 100.0]; % sigmas of gaussian kernel
t = 1.0;

img = phonecalc128; % try with phonecalc128, few128, nallo128
% im_out = gaussfft(im, 1.0);
psf1 = gaussfft(deltafcn(128,128),vars(1)); variance1 = variance(psf1);
psf2 = gaussfft(deltafcn(128,128),vars(2)); variance2 = variance(psf2);
psf3 = gaussfft(deltafcn(128,128),vars(3)); variance3 = variance(psf3);
psf4 = gaussfft(deltafcn(128,128),vars(4)); variance4 = variance(psf4);
psf5 = gaussfft(deltafcn(128,128),vars(5)); variance5 = variance(psf5);

figure()
subplot(1,5,1); showgrey(psf1); title('imp. response, t = 0.1');
subplot(1,5,2); showgrey(psf2); title('imp. response, t = 0.3');
subplot(1,5,3); showgrey(psf3); title('imp. response, t = 1.0');
subplot(1,5,4); showgrey(psf4); title('imp. response, t = 10.0');
subplot(1,5,5); showgrey(psf5); title('imp. response, t = 100.0');

img_out1 = gaussfft(img, t);
img_out2 = gaussfft(img, t*4);
img_out3 = gaussfft(img, t*16);
img_out4 = gaussfft(img, t*64);
img_out5 = gaussfft(img, t*256);

figure;
subplot(3,2,1); showgrey(img); title('Original image');
subplot(3,2,2); showgrey(img_out1); title('gauss filt im, t = 1.0');
subplot(3,2,3); showgrey(img_out2); title('gauss filt im, t = 4.0');
subplot(3,2,4); showgrey(img_out3); title('gauss filt im, t = 16.0');
subplot(3,2,5); showgrey(img_out4); title('gauss filt im, t = 64.0');
subplot(3,2,6); showgrey(img_out5); title('gauss filt im, t = 256.0');

%% test to build the function
img = phonecalc128;
vars = 100.0;
Im_in = fft2(img);  %compute the fourier transform
    
[dim1,dim2] = size(Im_in);
    
% [x,y] = meshgrid(0:dim1-1, 0:dim2-1);  %mesh for the Gaussian Kernel
x = -dim1/2:dim1/2 -1;
y = -dim2/2:dim2/2 -1;
    
% Gaussian Kernel
gauss = (1/2*pi*vars).*exp((-(x.^2  + y.^2)./(2*vars)));
% gauss = fftshift(gauss);

% gauss = (1/2*pi*var).*exp(-(x.^2)/(2*var)).*exp(-(y.^2)/(2*var));
% gauss = (1/2*pi*var).*exp((-((x-dim1/2).^2  + (y-dim2/2).^2)./(2*var)));
% gauss = fftshift(gauss);

%mesh(gauss);

%filtering
Im_out = Im_in.*gauss;
    
%output
im_out = ifft2(Im_out);

figure;
subplot(121); showgrey(img);
subplot(122); showgrey(im_out);

%% Part 3: Smoothing (Q17-Q18)
office = office256;
add = gaussnoise(office, 16); % add gaussian noise
sap = sapnoise(office, 0.1, 255); % add salt&pepper noise
ADD = fft2(add);
SAP = fft2(sap);
figure()
subplot(1,2,1); showfs(ADD); title("spectra SAP");
subplot(1,2,2); showfs(SAP); title("spectra SAP");


% 3 type of filters for each of noisy image
% for gaussian noisy, gaussan smooth, median and ideal LP
gauss_gauss = gaussfft(add, 0.9); %0.9
gauss_med = medfilt(add,3,3);
gauss_id = ideal(add,0.25);
% for sap noisy, gaussian smooth, median and ideal LP
sap_gauss = gaussfft(sap, 2); %2
sap_med = medfilt(sap,3,3);
sap_id = ideal(sap,0.2);
% plotting
figure()
subplot(3,3,1); showgrey(office); title('Original image');
subplot(3,3,2); showgrey(add); title('Image + gauss noise');
subplot(3,3,3); showgrey(sap); title('Image + sap noise');
subplot(3,3,4); showgrey(gauss_gauss); title('gauss-Gauss smooth (0.9)');
subplot(3,3,5); showgrey(gauss_med); title('gauss-Med filt (3x3)');
subplot(3,3,6); showgrey(gauss_id); title('gauss-Ideal LP filt (0.25)');
subplot(3,3,7); showgrey(sap_gauss); title('sap-Gauss smooth (2)');
subplot(3,3,8); showgrey(sap_med); title('sap-Med filt (3x3)');
subplot(3,3,9); showgrey(sap_id); title('sap-Ideal LP filt (0.2)');


%% Smoothing and subsampiing
img = phonecalc256;
%img = imread('wall.jpg');
%img = rgb2gray(img);
smoothimg = img;
smoothimg_ideal = img;
N = 5;
for i = 1 : 5
    if i > 1 % generate subsampled versions
        img = rawsubsample(img);
        smoothimg = gaussfft(smoothimg,2);
        smoothimg = rawsubsample(smoothimg);
        smoothimg_ideal = ideal(smoothimg_ideal,0.25);
        smoothimg_ideal = rawsubsample(smoothimg_ideal);
    end
    % subplot(3, N, i); imshow(img); title('Original subsampled'); % for image 'wall.jpg'
    subplot(3, N, i); showgrey(img); title('Original subsampled');
    % if i == 1 % for image 'wall.jpg'
        % subplot(3, N, i+N); imshow(smoothimg); title('Smoothed and subsampled');
        % subplot(3, N, i+2*N); imshow(smoothimg_ideal); title('Smoothed id and subsampled');
    % else % for image 'wall.jpg'
        subplot(3, N, i+N); showgrey(smoothimg); title('Smoothed (2) and subsampled ');
        subplot(3, N, i+2*N); showgrey(smoothimg_ideal); title('Smoothed id (0.25) and subsampled');
    % end

end


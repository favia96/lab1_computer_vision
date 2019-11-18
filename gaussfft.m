function [im_out, varMat] = gaussfft(im_in, var)
    Im_in = fft2(im_in);  %compute the fourier transform
    
    [dim1,dim2] = size(Im_in);
    
%     x = -dim1/2:dim1/2 -1;
%     y = -dim2/2:dim2/2 -1;
    
    [x,y] = meshgrid(0:dim1- 1, 0:dim1 -1);  %mesh for the Gaussian Kernel
    
    %Gaussian Kernel
    gauss = (1/2*pi*var).*exp((-((x-dim1/2).^2  + (y-dim2/2).^2)./(2*var)));
    gauss = fftshift(gauss);
%     gauss = (1/2*pi*var).*exp(-(x.^2)/(2*var)).*exp(-(y.^2)/(2*var));
%     figure;
%     mesh(gauss);
    
    %figure;
    %showgrey(ifft2(gauss));
    
    %filtering
    Im_out = Im_in.*gauss;
    
    %output
    im_out = ifft2(Im_out);
    varMat = variance(im_out);
end
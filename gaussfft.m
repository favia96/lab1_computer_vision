function im_out = gaussfft(im_in, var)
    Im_in = fft2(im_in);  % compute the fourier transform
    
    [dim1,dim2] = size(Im_in);
    
%   x = -dim1/2:dim1/2 -1;
%   y = -dim2/2:dim2/2 -1;
    
    [x,y] = meshgrid(0 : dim1- 1, 0 : dim2 -1);  %mesh for the Gaussian Kernel
    
    % Gaussian Kernel in freq
    gauss = (1/2*pi*var).*exp((-((x-dim1/2).^2  + (y-dim2/2).^2)./(2*var)));
    Gauss = fft2(gauss);
%   Gauss = exp((-((x-dim1/2).^2  + (y-dim2/2).^2).*(var/2)));
%   Gauss = fftshift(Gauss);

%   gauss = (1/2*pi*var).*exp(-(x.^2)/(2*var)).*exp(-(y.^2)/(2*var));
%   figure;
%   mesh(gauss);3.
    
    % filtering (multiplication in fourier)
    Im_out = Im_in .* Gauss;
    
    % output
    im_out = fftshift(ifft2(Im_out));
    % varMat = variance(gauss);

end
function[Fhat, F, Fabsmax] = fftwave(u, v, sz)
    if (nargin < 2)
        error('Requires at least two input arguments.')
    end
    if (nargin == 2)
        sz = 128;
    end

    Fhat = zeros(sz);
    Fhat(u, v) = 1;   %the only point equal to one in the frequency domain

    F = ifft2(Fhat);  %iverse Fourier tranform
%     F = fft2(Fhat);
    Fabsmax = max(abs(F(:)));   %max magnitude

    subplot(3, 2, 1);
    showgrey(Fhat);  %plot the image in the Fourier domain
    title(sprintf('Fhat: (u, v) = (%d, %d)', u, v))

    % these if and else move the origin of the image in the fourier domain
    % in the center of the image since when you compute direct transform
    % the center of the frequency space is exactly the center of the image,
    % in this section the new center is computed just for visualization
    % purpose since from the implementation pooint of view the frequency
    % image is shifted by using the instruction fftshift
    if (u <= sz/2)
        uc = u - 1;
    else
    uc = u - 1 - sz;
    end

    if (v <= sz/2)
    vc = v - 1;
    else
    vc = v - 1 - sz;
    end

%     wavelength = (2*pi)/sqrt((2*pi*u)^2 + (2*pi*v)^2); % Replace by correct expression
    wavelength = (2*pi)/sqrt(u^2 + v^2); % Replace by correct expression
    spectrum = sqrt(real(F).^2 + imag(F).^2);
    amplitude = max(spectrum(:)); % Replace by correct expression

    subplot(3, 2, 2);
    showgrey(fftshift(Fhat));
    title(sprintf('centered Fhat: (uc, vc) = (%d, %d)', uc, vc))

    subplot(3, 2, 3);
    showgrey(real(F), 64, -Fabsmax, Fabsmax);
    title('real(F)')

    subplot(3, 2, 4);
    showgrey(imag(F), 64, -Fabsmax, Fabsmax);
    title('imag(F)')

    subplot(3, 2, 5);
    showgrey(abs(F), 64, -Fabsmax, Fabsmax);
    title(sprintf('abs(F) (amplitude %f)', amplitude))

    subplot(3, 2, 6);
    showgrey(angle(F), 64, -pi, pi);
    title(sprintf('angle(F) (wavelength %f)', wavelength))
    
%     figure;
%     surf(real(F));
%     title('real part ');
%     
%     figure;
%     surf(imag(F));
%     title('imaginary part')
    
end
function[real_, Fhat, F, Fabsmax, amplitude, wavelength] = fftwave(u, v, sz)
    if (nargin < 2)
        error('Requires at least two input arguments.')
    end
    if (nargin == 2)
        sz = 128;
    end
    
    Fhat = zeros(sz); 
    Fhat(u, v) = 1; %image in Fourier domain, with an impulse in (u,v)

    F = ifft2(Fhat); %image in spatial domain, through ifft2, we have a sine
    Fabsmax = max(abs(F(:)));

    subplot(3, 2, 1);
    showgrey(Fhat);
    title(sprintf('Fhat: (u, v) = (%d, %d)', u, v))

    % shift the point (u,v) in the centered Fourier transform (0,0) in middle
    % of image
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

    real_ = real(F);
    wavelength = 2*pi / sqrt(u.^2 + v.^2) ; % Replace by correct expression
    spectrum = sqrt(real(F).^2 + imag(F).^2);
    amplitude = max(max(spectrum)); % Replace by correct expression
    
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
    
end
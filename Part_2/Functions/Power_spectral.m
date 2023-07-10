function y_fft_mag_squared = Power_spectral(y)

% Perform FFT on the signal
y_fft = fft(y)/length(y);

% Shift the spectrum to the center
y_fft_shifted = fftshift(y_fft);

% Calculate the magnitude squared
y_fft_mag_squared = abs(y_fft_shifted).^2;
end
clc;clear;
% Transmitter
bitstream_length = 100; % Length of the random bit stream
fc = 1e9; % Carrier frequency in Hz
Rb = 1;  % Bit rate (in this case, 1 bit per second)
ts=0.01;
 % Generate random bit stream
bitstream = randi([0, 1], 1, bitstream_length);
 % Line code the bit stream (Polar non return to zero)
linecoded_bits = 2 * bitstream - 1;
linecoded_bits = repelem(linecoded_bits, 1/ts);
 % Modulate the BPSK signal
t = linspace(0, bitstream_length/Rb, bitstream_length/ts); % Time vector
carrier = sqrt(2 * Rb) * cos(2 * pi * fc * t); % Carrier signal
modulated_signal = linecoded_bits .* carrier; % Modulated BPSK signal
% Spectral domain plot
figure(1);
y_fft_mag_squared = Power_spectral(linecoded_bits);
frequencies = (-fc/2:fc/10000:fc/2-fc/10000);
%frequencies = linspace(-1, 1, length(y_fft_shifted));
plot(frequencies, y_fft_mag_squared);
title('Spectrum of Line Coded BPSK Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude Squared');
grid on;
 % Plot the time domain of the modulated BPSK signal
 figure(2);
subplot(2,1,1);
plot(t, modulated_signal);
xlabel('Time');
ylabel('Amplitude');
title('Modulated BPSK Signal (Time Domain)');
 % Step 4: Plot the spectral domain using Fourier Transform
subplot(2,1,2);
plot(t,linecoded_bits);
xlabel('Time');
ylabel('Amplitude');
title('Line coded bits');
figure(3);
fs = 1/(t(2)-t(1)); % Sampling frequency
N = length(modulated_signal); % Length of the signal
freq = (-fc/2:fc/N:fc/2-fc/N); % Frequency vector
spectrum = Power_spectral(modulated_signal);
plot(freq, spectrum);
xlabel('Frequency');
ylabel('Magnitude');
title('Modulated BPSK Signal (Spectrum)');
% Receiver
received_signal = modulated_signal;
 % Demodulation (BPSK)
demodulated_signal = received_signal .* carrier;
 % Perform integration (LPF)
integration_factor = 0.01;
lpf_output = filter(integration_factor, [1, -integration_factor], received_signal);
% Set the threshold value
threshold = 0;

% Decision device before adding noise
for i = 1:length(lpf_output)
    if rem((i-1),1/ts)==0
        if lpf_output(i)>threshold
           received_bits((i+99)/(1/ts))=1;
        else
           received_bits((i+99)/(1/ts))=0;
        end
    end
end
 % Decision device
%received_bits = lpf_output > 0;
 % Count number of errors
num_errors = sum(bitstream ~= received_bits);

% Bit error rate (BER)
bit_error_rate = num_errors / bitstream_length;

disp(['Number of errors: ' num2str(num_errors)]);
disp(['Bit error rate: ' num2str(bit_error_rate)]);
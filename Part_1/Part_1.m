clc;clear;
%%%Set all variables

bitstream_length = 10000;       % Set the length of the bit stream
ts = 0.1;   % Set time samling
fs=1/ts;
t=0:ts:bitstream_length-ts;     % Set time vector
returnzero=repmat(repelem([1 0],1/(2*ts)),1,bitstream_length); % Set vector for RZ

%define eyediagram variables
n = 10;
period = 10;
offset = 6;

%define variables for the plot borders
%in case you want plot 10 bits set x_border to 10
%in case you want plot 10000 bits set x_border to 10000

x_border = 10;


%define frequency axis for power spectral
f = linspace(-fs/2,fs/2,bitstream_length*fs);

% generate 10 random values for sigma between 0 and 1.2
%sigma = rand(1, 10) * 1.2;
sigma = linspace(0,1.2,10);
%noise = sigma.*randn(1,length(bitstream_length));
% Generate a random bit stream
bitstream = randi([0,1], 1, bitstream_length);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified code to generate Unipolar NRZ/RZ with amplitude range from 0 to 1.2
ip = 1.2;
neg_ip = 0;
Unipolar=zeros(1,bitstream_length);

for i = 1:bitstream_length
    if bitstream(i) == 1
        Unipolar(i) = ip;
    else
        Unipolar(i) = neg_ip;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the Unipolar NRZ and plot it
UnipolarNRZ = repelem(Unipolar, 1/ts);
figure
plot(t, UnipolarNRZ);
axis([0 bitstream_length+1 -(ip+0.2) (ip+0.2)]);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Unipolar NRZ');
xlim([0 x_border])
ylim([-1.4 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Eye diagram
eyediagram(UnipolarNRZ,n,period,offset);
set(gcf, 'Name', 'Eye Diagram for Ideal Unipolar NRZ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot power spectral
y_fft_mag_squared= Power_spectral(UnipolarNRZ);

% Plot the spectrum
figure;
plot(f, y_fft_mag_squared);
xlabel('Normalized frequency (cycles/sample)');
ylabel('Power');
title('Power Spectral Density');

xlim([-5.00 5.00])
ylim([0.000000 0.000333])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Decesion device for Unipolar NRZ and check BER befor and after adding
%noise

% Initialize the received bitstream
received_bitstream = zeros(1, bitstream_length);

% Set the threshold value
threshold = 0.6;

% Decision device before adding noise
for i = 1:length(UnipolarNRZ)
    if rem((i-1),10)==0
        if UnipolarNRZ(i)>threshold
           received_bitstream((i+9)/10)=1;
        else
            received_bitstream((i+9)/10)=0;
        end
    end
end

% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) before adding noise
BER = num_errors/bitstream_length;

disp(['Bit Error Rate before adding noise for Unipolar_NRZ (BER): ' num2str(BER)]);

BER_Unipolar_NRZ =  Bit_Error_rate(UnipolarNRZ,sigma,bitstream,bitstream_length,t,threshold);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the Unipolar RZ and plot it
UnipolarRZ =UnipolarNRZ.*returnzero;
%t=0:ts:bitstream_length-ts;
figure;
plot(t, UnipolarRZ);
axis([0 bitstream_length+1 -(ip+0.2) (ip+0.2)]);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Unipolar RZ');
xlim([0 x_border])
ylim([-1.4 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Creating Eye diagram%%
eyediagram(UnipolarRZ,n,period,offset);
set(gcf, 'Name', 'Eye Diagram for Ideal Unipolar RZ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot power spectral
y_fft_mag_squared= Power_spectral(UnipolarRZ);

% Plot the spectrum
figure;
plot(f, y_fft_mag_squared);
xlabel('Normalized frequency (cycles/sample)');
ylabel('Power');
title('Power Spectral Density');

xlim([-5.00 5.00])
ylim([0.0000000 0.0000804])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Decesion device for Unipolar RZ and check BER befor and after adding
%noise

% Initialize the received bitstream
received_bitstream = zeros(1, bitstream_length);

% Set the threshold value
threshold = 0.6;

% Decision device before adding noise
for i = 1:length(UnipolarRZ)
    if rem((i-1),10)==0
        if UnipolarRZ(i)>threshold
           received_bitstream((i+9)/10)=1;
        else
            received_bitstream((i+9)/10)=0;
        end
    end
end

% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) before adding noise
BER = num_errors/bitstream_length;

disp(['Bit Error Rate before adding noise for Unipolar_RZ (BER): ' num2str(BER)]);

BER_Unipolar_RZ =  Bit_Error_rate(UnipolarRZ,sigma,bitstream,bitstream_length,t,threshold);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified code to generate polar NRZ/RZ with amplitude range from -1.2 to 1.2
ip = 1.2;
neg_ip = -1.2;
polar=zeros(1,bitstream_length);
for i = 1:bitstream_length
    if bitstream(i) == 1
        polar(i) = ip;
    else
        polar(i) = neg_ip;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the polar NRZ and plot it
polarNRZ = repelem(polar, 1/ts);
figure;
plot(t, polarNRZ);
axis([0 bitstream_length+1 -(ip+0.2) (ip+0.2)]);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('polar NRZ');
xlim([0 x_border])
ylim([-1.4 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Eye diagram
eyediagram(polarNRZ,n,period,offset);
set(gcf, 'Name', 'Eye Diagram for Ideal polar NRZ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot power spectral
y_fft_mag_squared= Power_spectral(polarNRZ);
% Plot the spectrum
figure;
plot(f, y_fft_mag_squared);
xlabel('Normalized frequency (cycles/sample)');
ylabel('Power');
title('Power Spectral Density');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Decesion device for Polar NRZ and check BER befor and after adding
%noise

% Initialize the received bitstream
received_bitstream = zeros(1, bitstream_length);

% Set the threshold value
threshold = 0;

% Decision device before adding noise
for i = 1:length(polarNRZ)
    if rem((i-1),10)==0
        if polarNRZ(i)>threshold
           received_bitstream((i+9)/10)=1;
        else
            received_bitstream((i+9)/10)=0;
        end
    end
end

% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) before adding noise
BER = num_errors/bitstream_length;

disp(['Bit Error Rate before adding noise for Polar_NRZ (BER): ' num2str(BER)]);

BER_Polar_NRZ =  Bit_Error_rate(polarNRZ,sigma,bitstream,bitstream_length,t,threshold);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the polar RZ and plot it
polarRZ =polarNRZ.*returnzero;
figure;
plot(t, polarRZ);
axis([0 bitstream_length+1 -(ip+0.2) (ip+0.2)]);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('polar RZ');
xlim([0 x_border])
ylim([-1.4 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Eye diagram
eyediagram(polarRZ,n,period,offset);
set(gcf, 'Name', 'Eye Diagram for Ideal polar RZ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot power spectral
y_fft_mag_squared = Power_spectral(polarRZ);

% Plot the spectrum
figure;
plot(f, y_fft_mag_squared);
xlabel('Normalized frequency (cycles/sample)');
ylabel('Power');
title('Power Spectral Density');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Decesion device for Polar RZ and check BER befor and after adding
%noise

% Initialize the received bitstream
received_bitstream = zeros(1, bitstream_length);

% Set the threshold value
threshold = 0;

% Decision device before adding noise to signal
for i = 1:length(polarRZ)
    if rem((i-1),10)==0
        if polarRZ(i)>threshold
           received_bitstream((i+9)/10)=1;
        else
            received_bitstream((i+9)/10)=0;
        end
    end
end

% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) before adding noise
BER = num_errors/bitstream_length;

disp(['Bit Error Rate before adding noise for Polar_RZ (BER): ' num2str(BER)]);



BER_Polar_RZ =  Bit_Error_rate(polarRZ,sigma,bitstream,bitstream_length,t,threshold);
% Modified code to generate Bipolar NRZ/RZ with amplitude range from -1.2 to 1.2
ip = 1.2;
neg_ip = 0;
Bipolar=zeros(1,bitstream_length);

for i = 1:bitstream_length
    if bitstream(i) == 1
        Bipolar(i) = ip;
        ip=-1*ip;
    else
        Bipolar(i) = neg_ip;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the bipolar NRZ and plot it
BipolarNRZ = repelem(Bipolar, 1/ts);
figure;
plot(t, BipolarNRZ);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Bipolar NRZ');
xlim([0 x_border])
ylim([-1.4 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Eye diagram
eyediagram(BipolarNRZ,n,period,offset);
set(gcf, 'Name', 'Eye Diagram for Ideal Bipolar NRZ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot power spectral
y_fft_mag_squared = Power_spectral(BipolarNRZ);

% Plot the spectrum
figure;
plot(f, y_fft_mag_squared);
xlabel('Normalized frequency (cycles/sample)');
ylabel('Power');
title('Power Spectral Density');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Decesion device for Bipolar NRZ and check BER befor and after adding
%noise

% Initialize the received bitstream
received_bitstream = zeros(1, bitstream_length);

% Set the threshold value
threshold_1 =0.6;
threshold_2 =-0.6;

% Decision device before adding noise to the signal
for i = 1:length(BipolarNRZ)
    if rem((i-1),10)==0
        if BipolarNRZ(i)>threshold_1 || BipolarNRZ(i)<threshold_2
           received_bitstream((i+9)/10)=1;
        else
            received_bitstream((i+9)/10)=0;
        end
    end
end

% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) before adding noise
BER = num_errors/bitstream_length;

disp(['Bit Error Rate before adding noise for Bipolar_NRZ (BER): ' num2str(BER)]);
BER_Bipolar_NRZ =  Bit_Error_rate(BipolarNRZ,sigma,bitstream,bitstream_length,t,threshold_1,threshold_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the Bipolar RZ and plot it
BipolarRZ =BipolarNRZ.*returnzero;
figure;
plot(t, BipolarRZ);
xlabel('Time (s)');
ylabel('Amplitude (V)')
title('Bipolar RZ');
xlim([0 x_border])
ylim([-1.4 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Eye diagram
eyediagram(BipolarRZ,n,period,offset);
set(gcf, 'Name', 'Eye Diagram for Ideal Bipolar RZ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot power spectral
y_fft_mag_squared = Power_spectral(BipolarRZ);

% Plot the spectrum
figure;
plot(f, y_fft_mag_squared);
xlabel('Normalized frequency (cycles/sample)');
ylabel('Power');
title('Power Spectral Density');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Decesion device for Bipolar RZ and check BER befor and after adding
%noise

% Initialize the received bitstream
received_bitstream = zeros(1, bitstream_length);

% Set the threshold value
threshold_1 =0.6;
threshold_2 =-0.6;

% Decision device before adding noise to the signal
for i = 1:length(BipolarRZ)
    if rem((i-1),10)==0
        if BipolarRZ(i)>threshold_1 || BipolarRZ(i)<threshold_2
           received_bitstream((i+9)/10)=1;
        else
            received_bitstream((i+9)/10)=0;
        end
    end
end

% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) before adding noise
BER = num_errors/bitstream_length;

disp(['Bit Error Rate before adding noise for Bipolar_RZ (BER): ' num2str(BER)]);

BER_Bipolar_RZ =  Bit_Error_rate(BipolarRZ,sigma,bitstream,bitstream_length,t,threshold_1,threshold_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified code to generate Manchester with amplitude range from -1.2 to 1.2
ip = 1.2;
neg_ip = -1.2;
Manchester=zeros(1,bitstream_length);
for i = 1:bitstream_length
    if bitstream(i) == 1
        Manchester(i) = ip;
    else
        Manchester(i) = neg_ip;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate the Manchester and plot it
Manchesterout=[];
one=repelem([1.2,-1.2],1/(2*ts));
zero=repelem([-1.2,1.2],1/(2*ts));
for i=1:bitstream_length
if Manchester(i)>0.6
Manchesterout=[Manchesterout one];
else
Manchesterout=[Manchesterout zero];
end
end
figure
plot(t, Manchesterout);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Manchester');
xlim([0 x_border])
ylim([-1.4 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating Eye diagram
eyediagram(Manchesterout,n,period,offset);
set(gcf, 'Name', 'Eye Diagram for Ideal Manchester');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot power spectral
y_fft_mag_squared = Power_spectral(Manchesterout);

% Plot the spectrum
figure;
plot(f, y_fft_mag_squared);
xlabel('Normalized frequency (cycles/sample)');
ylabel('Power');
title('Power Spectral Density');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Decesion device for Manchester and check BER befor and after adding
%noise

% Initialize the received bitstream
received_bitstream = zeros(1, bitstream_length);

% Set the threshold value
threshold = 0;

% Decision device before adding noise to signal
for i = 1:length(Manchesterout)
    if rem((i-1),10)==0
        if Manchesterout(i)>threshold
           received_bitstream((i+9)/10)=1;
        else
            received_bitstream((i+9)/10)=0;
        end
    end

end

% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) before adding noise
BER = num_errors/bitstream_length;

disp(['Bit Error Rate before adding noise for Manchester (BER): ' num2str(BER)]);

BER_Manchesterout =  Bit_Error_rate(Manchesterout,sigma,bitstream,bitstream_length,t,threshold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting BERs with sigma

% with smoothing curves

smoothed_BER_Manchesterout = movmean(BER_Manchesterout,3);
smoothed_BER_Bipolar_RZ = movmean(BER_Bipolar_RZ,3);
smoothed_BER_Bipolar_NRZ = movmean(BER_Bipolar_NRZ,3);
smoothed_BER_Polar_RZ = movmean(BER_Polar_RZ,3);
smoothed_BER_Polar_NRZ = movmean(BER_Polar_NRZ,3);
smoothed_BER_Unipolar_RZ = movmean(BER_Unipolar_RZ,3);
smoothed_BER_Unipolar_NRZ = movmean(BER_Unipolar_NRZ,3);
figure
semilogy(sigma,smoothed_BER_Manchesterout)
hold on
semilogy(sigma,smoothed_BER_Bipolar_RZ)
semilogy(sigma,smoothed_BER_Bipolar_NRZ)
semilogy(sigma,smoothed_BER_Polar_RZ)
semilogy(sigma,smoothed_BER_Polar_NRZ)
semilogy(sigma,smoothed_BER_Unipolar_RZ)
semilogy(sigma,smoothed_BER_Unipolar_NRZ)
legend('Manchester Out', 'Bipolar RZ', 'Bipolar NRZ', 'Polar RZ', 'Polar NRZ', 'Unipolar RZ', 'Unipolar NRZ','Location','southeast')

grid on
title('BER versus Sigma')
xlabel('Sigma')
ylabel('BER')


% without smoothing the curve
%{
semilogy(sigma,BER_Manchesterout)
hold on
semilogy(sigma,BER_Bipolar_RZ)
semilogy(sigma,BER_Bipolar_NRZ)
semilogy(sigma,BER_Polar_RZ)
semilogy(sigma,BER_Polar_NRZ)
semilogy(sigma,BER_Unipolar_RZ)
semilogy(sigma,BER_Unipolar_NRZ)
legend('Manchester Out', 'Bipolar RZ', 'Bipolar NRZ', 'Polar RZ', 'Polar NRZ', 'Unipolar RZ', 'Unipolar NRZ','Location','southeast')
%}
%Bonus part
%  adding noise  and counting the number of errors for the bits
for i=1:10
threshold_1=0.9;
threshold_2=-0.9;
noise = sigma(i)*randn(1,length(t));
BipolarRZ_noise = noise + BipolarRZ ;

%taking the average value for the noise ber 10 bits
BipolarRZ_noise_average=zeros(1,length(BipolarRZ_noise));
for q=1:bitstream_length
BipolarRZ_noise_average(10*(q-1)+1:10*(q-1)+5)=mean(BipolarRZ_noise(10*(q-1)+1:10*(q-1)+5));
end
 errors = 0;
last_one=-1.2;
 for j = 1:length(BipolarRZ_noise_average)-11
  if rem((j-1),10)==0
        if BipolarRZ_noise_average(j) > threshold_1 && last_one== 1.2
        errors = errors + 1;
        last_one=1.2;
        elseif BipolarRZ_noise_average(j) >threshold_1  && last_one== -1.2
        last_one=1.2;
        elseif BipolarRZ_noise_average(j) < threshold_2 && last_one== -1.2
        errors = errors + 1;
        last_one=-1.2;
        elseif BipolarRZ_noise_average(j) < threshold_2 && last_one== 1.2
        last_one=-1.2;

       end
  end
end
disp(['Number of Errors out of Sigma' num2str(i)]);
disp([' = ' num2str(errors)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the Bipolar RZ with the noise after the detection
figure
    subplot(2,2,1)
    plot(t, BipolarRZ_noise_average);
    xlabel('Time (s)');
    ylabel('Amplitude (V)')
    title('Bipolar RZ with noise ');
    xlim([0 10]);
    ylim([-3 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % plot the Bipolar RZ without noise after the detection
    subplot(2,2,2)
    plot(t, BipolarRZ);
    xlabel('Time (s)');
    ylabel('Amplitude (V)')
    title('Bipolar RZ without noise ');
    xlim([0 10]);
    ylim([-3 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % plot the Streamof bits at start
 output_signal = repelem(bitstream, 1/ts);
    subplot(2,2,3)
    plot(t, output_signal);
    xlabel('Time (s)');
    ylabel('Amplitude (V)')
    title('The stream of bits at start ');
    xlim([0 10]);
    ylim([-3 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision device after adding noise to the signal
threshold_1 =0.9;
threshold_2 =-0.9;
for ii = 1:length(BipolarRZ_noise_average)
    if rem((ii-1),10)==0
        if BipolarRZ_noise_average(ii)>threshold_1 ||  BipolarRZ_noise_average(ii)<threshold_2
           received_bitstream((ii+9)/10)=1;
        else
            received_bitstream((ii+9)/10)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Decision Decive to detect the BIpolar RZ with noise
 received_signal = repelem(received_bitstream, 1/ts);
    subplot(2,2,4)
    plot(t,received_signal);
    xlabel('Time (s)');
    ylabel('Amplitude (V)')
    title('The stream of bits recived ');
    xlim([0 10]);
    ylim([-3 3]);
end
function BER_Unipolar_NRZ =  Bit_Error_rate(x,sigma,bitstream,bitstream_length,t,threshold_1,threshold_2)
    
if nargin == 6
for i = 1:length(sigma)
    noise = sigma(i)*randn(1,length(t));
    x_noise = noise + x ;
    for j = 1:length(x_noise)
    if rem((j-1),10)==0
        if x_noise(j)>threshold_1 
           received_bitstream((j+9)/10)=1;
        else
           received_bitstream((j+9)/10)=0;
        end
    end
    end
% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) after adding noise
BER = num_errors/bitstream_length;
BER_Unipolar_NRZ(i)=BER;
disp(['Bit Error Rate after adding noise for the signal (BER): ' num2str(BER)]);
end
elseif nargin ==7
    for i = 1:length(sigma)
    noise = sigma(i)*randn(1,length(t));
    x_noise = noise + x ;
    for j = 1:length(x_noise)
    if rem((j-1),10)==0
        if x_noise(j)>threshold_1 || x_noise(j)<threshold_2
           received_bitstream((j+9)/10)=1;
        else
           received_bitstream((j+9)/10)=0;
        end
    end
    end
% Compare with the transmitted bitstream and count errors
num_errors = sum(bitstream ~= received_bitstream);

% Calculate bit error rate (BER) after adding noise
BER = num_errors/bitstream_length;
BER_Unipolar_NRZ(i)=BER;
disp(['Bit Error Rate after adding noise for the signal (BER): ' num2str(BER)]);
    end
end
end


    

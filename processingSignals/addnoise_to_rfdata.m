function [ noise_data ] = addnoise_to_rfdata( rfdata,  snr)
[num_sample,num_receiver,num_transmitter] = size(rfdata);
noise_data = zeros(num_sample,num_receiver,num_transmitter);
for ii = 1:num_transmitter
    noise_data(:,:,ii)  = addNoise(rfdata(:,:,ii),snr);
end
end


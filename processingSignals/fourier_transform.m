function fourier_transform (wave, sampling_freq)
fs = sampling_freq;
y = fft(wave);
n = length(wave);
f = (0:n-1) * (fs/n);
power = abs(y).^2/n;

figure;
plot(f,power)
xlabel('Frequency')
ylabel('Power')

y0 = fftshift(y);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(y0).^2/n;    % 0-centered power

figure;
plot(f0,power0)
xlabel('Frequency')
ylabel('Power')

end
function [ noisy, time ] = EIT_add_colored_noise( sig, PSD )
%EIT_ADD_COLORED_NOISE Add coloured noise according to the PSD of the
%signal
%   (C) 2015/02/18 Sebastien Martin

debug = false;
% 1st part from awgn function
sigPower = sum(abs(sig(:)).^2)/length(sig(:));
sigPower = 10*log10(sigPower);
noisePower = sigPower-max(abs(PSD));
% This is from wgn function, called from awgn
noisePower = 10^(noisePower/10);
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/336939
[noise, time] = get_noise(PSD, false); % No need to debug
noise = sqrt(noisePower) * noise(1:length(PSD));

noisy = sig + noise;

if debug;
    figure; plot([sig', noisy']); legend('Initial signal', 'Noisy Signal');
    fft_sig = fft(sig);
    fft_sig = fftshift(abs( fft_sig )) / size(fft_sig,1);
    fft_noisy = fft(noisy);
    fft_noisy = fftshift(abs( fft_noisy )) / size(fft_noisy,1);
    
    freq_max = 1000; % Sample at 2MHz, maximal frequency measurable (without aliasing and stuffs)
    freqs = linspace(0, freq_max, size(fft_noisy,2)/2);
    fft_plot = [fft_sig; fft_noisy]';
    figure; plot(freqs, fft_plot(size(fft_plot,1)/2+1:end,:));
    xlabel('Frequency (kHz)');
    ylabel('Amplitude spectrum');
    legend('Initial signal (FFT)', 'Noisy Signal (FFT)');
    title('FFT of simulated signal');
end
end

function [ waveTs, time ] = get_noise( PSD, debug )
% Create timeseries from a given double-sided PSD
%
% Create perfect white noise by generating random phase in the frequency
% domain, and then multiply it by your PSD. Then transfer it to the time
% domain...Take note that MATLAB expects the conjugate symmetric part of
% the spectrum to be the second half of the array.
%
% INPUTS
% PSD - Double-sided power spectral density in the expected MATLAB
% order e.g. [wavePSD_positiveFreq wavePSD_negativeFreq]
%
% fs - Sample rate of the output timeseries [samples/sec]
%
% T - Length of the output timeseries [seconds]
%
% OUTPUTS
% waveTs - Generated timeseries
%
% time - Time vector corresponding to timeseries generated
%
%
% Copyright Mike Rudolph - 2014

T = (length(PSD)-1)*2; % Length of your timeseries [seconds]
fs = 2000000; % Sample rate [samples/sec]
T = T/fs;

dt = 1/fs; % Delta time [sec]
df = 1/T; % Delta frequency [Hz]
N = round(fs*T);

frequency_DoubleSided = (0:(N-1))*df; % Double-Sided Frequency vector
time = (0:N-1)*dt; % Time vector


%% Your PSD, should be the same length as noise
% Make it double-sided
wavePSD_positiveFreq = PSD;
wavePSD_negativeFreq = flip(conj(wavePSD_positiveFreq(2:end-1)));
wavePSD_DoubleSided = [wavePSD_positiveFreq wavePSD_negativeFreq];

%% Generate 'Perfect' white noise in the frequency domain
if rem(N,2) == 1
    N = N-1;
end
randnums = rand(1, N/2-1).*2*pi; % Random phase between 0 and 2pi
randvalues = exp(1i*randnums); % This is your white noise
linspecPositiveFreq = [1 randvalues 1]; % Positive Frequencies
linspecNegativeFreq = flip(conj(randvalues)); % Negative Frequencies

% Need this order for IFFT in MATLAB:
noisePSD_DoubleSided = [linspecPositiveFreq linspecNegativeFreq];


%% Multiply noise by PSD (both double-sided) in frequency-domain
totalWavePSD = wavePSD_DoubleSided.*noisePSD_DoubleSided;

% Convert double-sided PSD to time domain via IFFT
waveTs = real(ifft(totalWavePSD)); % Should be all real anyway


%% Plot it up, yo.
if debug
    figure;
    subplot(211)
    plot(frequency_DoubleSided, abs(totalWavePSD))
    title('Double-Sided PSD * White Noise Magnitude')
    xlabel('Frequency (Hz)')
    ylabel('Spectral Density Magnitude [m^2/Hz]')
    
    subplot(212)
    plot(frequency_DoubleSided, angle(totalWavePSD))
    title('Double-Sided PSD * White Noise Phase')
    xlabel('Frequency (Hz)')
    ylabel('Phase [rad]')
    
    figure;
    plot(time, waveTs)
    title('Generated Timeseries From PSD')
    xlabel('Time (Sec)')
    ylabel('Amplitude [m]')
end

end
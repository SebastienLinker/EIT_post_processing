function [ SNR ] = EIT_meas_est_noise( raw_sig )
%EIT_MEAS_EST_NOISE Estimate noise in EIT measurements
%   Assumes a white gaussian noise and estimae the SNR of the signal
%	Input: Noisy signals
%	Output the estimated SNR, in dB, for each measurement (same electrode
%	for injection, same electrode for measurement)
%	see http://dsp.stackexchange.com/questions/4889/how-do-i-calculate-snr-of-noisy-signal

% Note: Requires Symbolic Math Toolbox

debug = false;

% Divide data
n_inj = size(raw_sig,1)/60; % 3*20samples
n_elec = size(raw_sig,2);
raw_meas = reshape(raw_sig,60,n_inj*n_elec); % Individual sine waves

% Assume the expected signal is a perfect sine wave
amps = max(raw_meas);
% Estimate highest peak and shift the sine wave
expect = zeros(60, n_elec*n_inj);
sin_pt = [0:1:59]/60*(2*pi)*3; % Perfect sine wave
for k=1:1:n_elec*n_inj
    % Position of maximal amplitude
    sin_sh = asin(raw_meas(1,k))/abs(amps(k)); % Shift
    if ((raw_meas(1,k)<0) && (raw_meas(1,k)>raw_meas(2,k)))
        sin_sh=-sin_sh;
    end
    l=1;
    while (sign(raw_meas(l,k)) ~= sign(raw_meas(l+1,k))); l=l+1; end
    sin_dir = heaviside(abs(raw_meas(l+1,k))-abs(raw_meas(l,k)))*pi; % Direction
    if raw_meas(l,k)>0; sin_dir = sin_dir + pi;	end
    
    expect(:,k) = (abs(amps(k)) * sin( sin_pt + sin_sh + sin_dir) )'; % Expected sine wave
    if debug
        figure; plot([expect(:,k), raw_meas(:,k)]);
        disp(['Expected 1st value: ',num2str(raw_meas(1,k)),' 1st value: ',num2str(expect(1,k))]);
        disp(['Difference_first two signals: ',num2str(abs(raw_meas(2,k))-abs(raw_meas(1,k)))]);
    end
end

res_n = raw_meas - expect; % Residual noise
m_raw_meas = mean(raw_meas.^2);
m_res_n = mean(res_n.^2);

SNR = zeros(n_inj*n_elec,1);
for k=1:1:n_inj*n_elec
    SNR(k) = m_raw_meas(:,k) / m_res_n(:,k);
end
SNR = 10*log10(SNR);

if debug;
    figure; plot(SNR); title('Estimated SNR on noisy data');
    disp(['Estimated noise average: ',num2str(mean(SNR)),' std: ',num2str(std(SNR))]);
end

end
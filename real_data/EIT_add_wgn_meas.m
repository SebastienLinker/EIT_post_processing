function [ EIT_meas ] = EIT_add_wgn_meas( EIT_meas, noise, filt, n_elecs)
%EIT_ADD_WGN_MEAS Add Noise to measurements
%   The noise on measurements is assumed to be a White Gaussian Noise,
%   filtered.
%	The filter is assumed to be a bandpass butterworth filter, order 10,
%	and is similar to the filter used in our hardware measurements
%	(phantom). The hardware injects AC current at 100kHz and measures at
%	2MHz. For each injection, 3 sine waves are acquired, duplicated,
%	filtered, and the maximal amplitude of each sine wave is detected. All
%	the measurements are done at these 3 points, and averaged.
%	@input
%	EIT_meas: The set of measurements
%	noise: The noise model
%	filt: boolean, filter data or not (default: true)

% Note: Requires Communications toolbox

debug = false;

if nargin < 3; filt=true; end
if nargin < 4; n_elecs=8; end

%% Go back to 3 sine waves
n_meas_p_inj = length(EIT_meas)/n_elecs;
max_amp = EIT_meas;
sin_pt = [0:1:59]/60*(2*pi)*3; % Perfect sine wave
s_waves = (max_amp * sin(sin_pt))';
if debug
    figure; plot(s_waves); title('Back to sine waves');
    figure; plot(noise); ylabel('SNR (dB)'); title('Noise model');
end

%% Add noise
tot_len = length(EIT_meas);
noisy_sin = zeros(60, tot_len);
for k=1:1:tot_len
    if isscalar(noise(k,:))
        noisy_sin(:,k) = awgn(s_waves(:,k), noise(k), 'measured');
    else
        noisy_sin(:,k) = EIT_add_colored_noise(s_waves(:,k)', noise(k,:))';
    end
    % 	if debug; figure; plot(noisy_sin(:,k)); title('Measurements with noise (Simulation)'); end
end

%% Filter the data
if filt
    filtered = filter_EIT_meas(noisy_sin, n_elecs);
else
    filtered = noisy_sin;
end
if debug
    dbg_filt = reshape(filtered, size(noisy_sin));
    for k=1:1:length(EIT_meas)
        figure;
        plot([noisy_sin(:,k),dbg_filt(:,k)]);
        title('Filtered vs unfiltered');
    end
end

%% Re estimate the noise (debug purpose)
if debug
    n_inj = tot_len/n_elecs;
    est_noisy = EIT_meas_est_noise( reshape(noisy_sin, n_inj*60, n_elecs) );
    est_filt = EIT_meas_est_noise( reshape(filtered, n_inj*60, n_elecs) );
    figure; plot([est_noisy, est_filt, noise]);
    xlabel('Measurements'); ylabel('noise (dB)');
    legend('No filter','Filtered','Noise model');
    title('Noise estimate (from simulation)');
end

%% Go back to the EIT measurement (EIDORS format)
EIT_meas = detect_Amp(filtered);
% Error correction
to_switch = (sign(EIT_meas)~=sign(max_amp));
EIT_meas(to_switch) = -EIT_meas(to_switch);

if debug
    figure; plot([max_amp, EIT_meas]); legend('No noise', 'Noise');
end

end
function [ filtered ] = filter_EIT_meas( unfiltered, n_elecs )
%FILTER_EIT_MEAS Simulate a filter on EIT measurement
%	For each injection, the signal is acquired over 3 periods
%   Sine wave: 100 kHz
%	Measurement: 2 MHZ
%	Order 10, cut-off frequencies: 95kHz, 105kHz

% Note: Requires Signal processing toolbox

debug = false;

if nargin==1
    n_elecs=8;
end

sine_fr = 100000;
fs = 2000000;
wind_fr = 5000; % Window
Att_pb = 0.025; %Attenuation in pass band
Att_sb = 90; % Attenuation in stop band
order = 10; % Butter function multiplies by 2 for passband filters

if debug; figure; plot(unfiltered); end

Hd = eidors_cache(@get_filter, {sine_fr, fs, wind_fr, Att_pb, Att_sb}, 'phantom_filter');

if debug; fvtool(Hd,'Fs',fs); end

% Divide data
n_inj = size(unfiltered,1)/60; % 3*20samples
n_meas = size(unfiltered,2);
to_filt = reshape(unfiltered,60,n_inj*n_meas);
to_filt = repmat(to_filt,100,1);

% Filter
filt = filter(Hd,to_filt);
if debug
    figure; plot(to_filt); title('To filter');
    figure; plot(filt); title('Filtered');
end

% Choose the last 3 periods
filtered = filt([end-(59*2):1:end-59],:);

% Re-order everything
filtered = reshape(filtered, size(filtered,2)*60/n_meas, n_meas);

if debug; figure; plot(filtered); end

end

function Hd = get_filter(sine_fr, fs, wind_fr, Att_pb, Att_sb)
% Create filter
Cf1 = sine_fr-wind_fr;
Cf2 = sine_fr+wind_fr;
% [b, a] = butter( order, [Cf1 Cf2]/fs, 'bandpass');

Sf1 = Cf1-wind_fr;
Sf2 = Cf2+wind_fr;
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', Sf1, Cf1, Cf2, Sf2, Att_sb,Att_pb,Att_sb, fs);
Hd = design(d,'butter');
end
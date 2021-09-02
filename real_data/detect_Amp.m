function [ EIT_meas, var_inj ] = detect_Amp( s_wave )
%DETECT_AMP For each sine wave, detect the highest amplitude
%   Detailed explanation goes here

debug = false;
if debug; figure; plot(s_wave); title('Sine wave before detect amp'); end

% Divide data
n_meas = size(s_wave,1)/60; % 3*20samples
n_elec = size(s_wave,2);
s_meas = reshape(s_wave,60,n_meas*n_elec);
per_inj = reshape(s_meas,60,n_elec,n_meas);

if debug; for k=1:n_meas; figure; plot(per_inj(:,:,k)); title(['Injection ',int2str(k)]);end;  end

% Detect maximal or minimal amplitude
I = zeros(n_meas,60);
amp_inj = zeros(n_elec,n_meas);
var_inj = zeros(n_elec,n_meas);
[~, pk] = max(abs(per_inj(:)));
dir = 'descend'; if per_inj(pk)<0; dir='ascend'; end
for k = 1:1:n_meas
    [~, I(k,:)] = sort(per_inj(:,k,k),dir);
    % Do the measurements when the injection is at the maximum (peak of sine wave)
    % 	I(k,1:3) = [15 35 55];
    amp_inj(:,k) = mean(per_inj(I(k,1:3),:,k));
    var_inj(:,k) = var(per_inj(I(k,1:3),:,k));
end

EIT_meas = reshape(amp_inj,n_meas*n_elec,1);

if debug
    figure; plot(amp_inj'); title('Highest amplitude detected');
    figure; plot(EIT_meas); title('EIT full set of measurements');
end

end
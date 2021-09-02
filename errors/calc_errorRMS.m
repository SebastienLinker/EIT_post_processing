function [ RMS ] = calc_errorRMS( target, reconstructed, varargin )
%CALC_ERRORRMS Compute the error, RMS method
%   Compute the RMS error
%	Output the RMS error
%	Please refer to function calc_error_MSE for further documentation
%	@Input :
%		target
%		reconstructed
%		varargin: additional arguments, can be one of the following strings
%			area: compute the size of the elements first, useful if the
%				model has large and small elements
%			volume: same as area, for 3D models
%			region: Equation specifying the error region
%			list: List of elements on which we should calculate the error
%			normalize: In case the models are not normalized
%			...
%	@Output :
%		RMS: The resulting RMS error
%
%	See also calc_Error_MSE

if any(size(target)~=size(reconstructed))
    target = repmat(target, size(reconstructed)./size(target));
end

totlen = size(reconstructed);
RMS = zeros(totlen);

for k=1:1:totlen(1)
    for l=1:1:totlen(2)
        [~, MSE] = calc_error_MSE( target(k,l), reconstructed(k,l), varargin{:} );
        RMS(k,l) = sqrt(MSE);
    end
end

end
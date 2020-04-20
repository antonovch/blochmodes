function [ld0, Q] = BICModel(varargin)
% function	[ld0, Q] = BICModel(kx,lambda,data,h,bw_cutoff,max_p)
%	for homogeneous media, symmetric
% function	[ld0, Q] = BICModel(kx,lambda,data,data2,h,bw_cutoff,max_p)
%	for substrate, different top and bottom
%
% Applies either 2 Bloch wave model with 2 BW as the dominant one,
% or 3 BW model with either BW 2 or 3 as the dominant
%
% INPUT ARGS:
%   - kx: in-plane wavevector (unnormalized, scalar/array)
%   - lambda: wavelength values (scalar/array)
%   - data: struct array with fields being reflection/transmission coeffs and n_eff
%       each struct corresponds to one value of kx
%   - data2: structure for the second interface in the broken vertical symmetry case
%   - h: slab thinkness (scalar/array)
%   - bw_cutoff(1): selects both the model and lambda cutoff (this often
%   makes no difference). The number was supposed to specify the BW cutoff 
%   up tp which we're calculating
%       4 -- 3 BW model with the 3rd as the dominant one
%       4i -- 3 BW with the 2nd as the dominant one
%       3 -- 2 BW with the 2nd as the dominant one
%       3i -- 2 BW, 1 and 3, without the second one
%       bw_cutoff can be a vector of two elements to indicate start/stop cutoffs
%       the stop cutoff (second element in bw_cutoff) allows to narrow the band
%       for interpolation
%   - max_p: maximum phase-matching order
%
% OUTPUT ARGS:
%   - ld0: PhC modes eigenwavelengths (cell array of length max_p+1)
%   - Q: Q-factors of the PhC modes (cell array of length max_p+1)
%
% Not implemented: effective transmission for 3 BW model

if length(varargin)==6
	[kx,lambda,data,h,bw_cutoff,max_p] = deal(varargin{:});
elseif length(varargin)==7
	[kx,lambda,data,data2,h,bw_cutoff,max_p] = deal(varargin{:});
else
	error('Wrong number of arguments')
end

reverseStr = ''; % to display the progress
switch bw_cutoff(1)
	case {4,4i,3i} % including the region with 3 BWs, not just 2
		ld_start = cal_cutoff(data,4);
	case 3 % including only the region with 2 BWs 
		ld_start = cal_cutoff(data,3);
	otherwise
		error('Wrong bw_cutoff')
end

if length(bw_cutoff)==2
	switch bw_cutoff(2)
		case 3 % only the region with 3 BWs, not just 2
			ld_stop = cal_cutoff(data,3);
		case 2 % also region with 2 BWs 
			ld_stop = cal_cutoff(data,2);
		otherwise
			error('Wrong bw_cutoff')
	end
else
	ld_stop = ones(1,length(kx))*length(lambda);
end

[ld0, Q, Reff0, dPhiT0] = deal(cell(1,max_p+1));
[ld0{:}, Q{:}, Reff0{:}, dPhiT0{:}] = deal(zeros(length(kx),length(h)));

dld = 1e-4; % step for finite difference derivative

for j=1:length(kx)
	if length(varargin)==6
		[reff, phiT, teff] = calcEffectiveCoeffsSym(lambda,data(j),h,bw_cutoff(1),[ld_start(j) ld_stop(j)]);
		Reff=abs(reff.^2);
	elseif length(varargin)==7
		[reff, reff2, phiT, teff, teff2] = calcEffectiveCoeffsAsym(lambda,data(j),data2(j),h,bw_cutoff(1),[ld_start(j) ld_stop(j)]);
		Reff=abs(reff.*reff2);
    end
    ld = lambda(ld_start(j):ld_stop(j));
	for ii=1:length(h)
	   tmpPhiT=phiT(:,ii)/pi; 
       % rarely, non-unique, values, notably, zeros, but also nans, can occur,
       % which cause interp1 to fail, so we remove them
	   [~, ia, ~] = unique(tmpPhiT,'first'); 
       I = zeros(size(tmpPhiT));
       I(ia) = 1;
       I = logical(I) & ~isnan(tmpPhiT);
	   for ord=0:max_p
           if max(tmpPhiT)>ord && min(tmpPhiT)<ord
               ld0{ord+1}(j,ii)=interp1(tmpPhiT(I),ld(I),ord,'spline');
               Reff0{ord+1}(j,ii)=interp1(ld(I), Reff(I,ii), ld0{ord+1}(j,ii),'spline');
               dPhiT0{ord+1}(j,ii)=(interp1(ld(I), phiT(I,ii), ld0{ord+1}(j,ii)+dld,'spline')-...
               	interp1(ld(I), phiT(I,ii), ld0{ord+1}(j,ii)-dld,'spline'))/(2*dld);
           end
	   end
	end
	for ord=0:max_p
	    Q{ord+1}(j,:)=-ld0{ord+1}(j,:)./(1-Reff0{ord+1}(j,:)).*dPhiT0{ord+1}(j,:);
	end

	% Display the progress
   msg = sprintf('kx values done: %g/%g\n', j, length(kx));
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

end

function n_eff_cutoff = cal_cutoff(data, n)

	n_eff_cutoff = ones(1,length(data));
	if size(data(1).n_eff,1) < n
		return
	end
	for j = 1:length(data)
		if ~isempty(find(data(j).n_eff(n,:)==0,1))
		    n_eff_cutoff(j)=find(data(j).n_eff(n,:)==0,1)-1;
		else
		    n_eff_cutoff(j)=length(data(j).n_eff(n,:));
		end
	end
end
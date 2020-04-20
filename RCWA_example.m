
%%%% Parameters
nfourier = 30;  % Fourier harmonics run -nfourier:nfourier
ff = 0.6; % filling factor in dielectric material
period = 1;    % period in µm, equals to unity, can be viewed as evertyhing is normatized to p
nd = 3.5; % refractive index of the dielectric material
n_air = 1.0; % refractive index in the grating slits and above the grating

%%%% Sweep variables
% all lengths has to be of the same units, here it's microns
lambda = 2.4; 
kx = 0.23... % x-axis units in the paper
    *2*pi/period... % normalization from the paper
    /(2*pi./lambda); % normalization to k0 for res1 function

%%%% Reticolo definitions
parm = res0(1); % default parameters for TE polarization, 
                % TM : parm=res0(-1), see Reticolo documentation 
parm.res1.trace = 0; % plot the layers index profile
parm.res1.ftemp = 0; % do not create temp files (default = 1), delete with retio call 

textures=cell(1,2); % see Reticolo documentation for details on texture definition
textures{1} = n_air; % uniform above 
textures{2} = {[-ff/2, ff/2], [n_air, nd]}; % grating region

% a - stores layer information for S-matrix calculation
% n_eff - cell array with modes effective indices for each layer. Number of
% modes = 2*nfourier+1
[a, n_eff] = res1(lambda, period, textures, nfourier, kx, parm);

%%%% S-Matrix above
[s_top, top, k_par_top] = retb(a{2}, a{1}{1}, 1e-6);
% truncation of non-propagative modes
f = find(abs(kx-k_par_top(:,1))<100*eps);
s_top = rettronc(s_top, top(f,1), top(f,1), 1);

%%%% S-matrix below and in the grating (for the half-problem)
[s_bottom, bottom, k_par_bottom] = retb(a{2}, a{1}{2}, -1e-6);
% truncation of non-propagative modes
f = find(abs(real(bottom(:,2))) == 0);
[~, I] = sort(abs(imag(bottom(f,2))), 1, 'descend');
s_bottom = rettronc(s_bottom, bottom(f(I),1), bottom(f(I),1), -1);

s_inter = retss(s_top, s_bottom); 
S = s_inter{1};

% retio % erase temporary files

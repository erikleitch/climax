addpath matlab
addpath /home/szadaq/sza_analysis/matutils
addpath /home/szadaq/sza_analysis/matutils/interf
addpath /home/szadaq/sza_analysis/constants

'Calculating beams...'

ad=calc_ad(0.5, 256);
ad.nu       = sza_bandfreq(1:2);

ad.Tnorm_uK = sqrt(500); % Normalize to bandpower of 500 uK^2
ad.lnorm = 4000;         % At l = 4000
ad.type = 'pow';
ad.powInd = -2;          % Make the spectrum of C_l -2, so that the bandpower 
                         % (=(l*(l+1)*C_l ) is approximately flat

ad.beam = ones(256,256,2);

dx   = (ad.Field_size_deg/ad.N_pix)/180*pi;

%du   = 1.0/(ad.N_pix * dx) / 2;
%umax = 1.0/(2 * dx) / 2;

du   = 1.0/(ad.N_pix * dx);
umax = 1.0/(2 * dx);

us = -umax:du:umax;
vs = -umax:du:umax;

uss(:,:,1) = repmat(us,257,1);
vss(:,:,1) = repmat(vs,257,1)';

uss(:,:,2) = repmat(us,257,1);
vss(:,:,2) = repmat(vs,257,1)';

fs = ones(257,257,2);
fs(:,:,1) = ad.nu(1);
fs(:,:,2) = ad.nu(2);

vis = complex(uss,vss);
var = complex(fs,fs);

szaMatSimVis(ad, uss, vss, vis);

lss = sqrt(uss .^2 + vss .^2);

ls = lss(:,:,1);
rs = real(vis(:,:,1));

dl = 500;
nbin = 10000/dl;

vvs = zeros(1,nbin);
lvs = zeros(1,nbin);

for i=1:nbin
  lmin = (i-1) * dl;
  lmax = (i) * dl;

  inds  = find(ls > lmin & ls <= lmax);
  rvals = rs(inds);
  lvals = ls(inds);

  vvs(i) = rms(rvals(:)).^2;  
  lvs(i) = mean(lvals(:));
end


%------------------------------------------------------------
% Main method of this file.  Read in the miriad file, rearrange the
% data, and write it out.  Use like:
%
%     writeuvf(mirfile, srcname, uvffile, calcWt)
%
% Where:
%
%   mirfile = input miriad file
%   srcname = source to extract
%   uvffile = output uvf file
%   calcWt  = calculate weights from Tsys (1) or estimate from data (0)
%
%------------------------------------------------------------

function writeuvf(mirfile, srcname, uvffile, calcWt)

  d = climaxMatReadMiriadEml(mirfile);

  [g coord] = getSrc(d, srcname);

  [nAnt antmap revmap bases] = getBaselineMapping(g)

  nBase = (nAnt*(nAnt-1))/2;

  [data uvw jd freq bw] = getData(d, nAnt, nBase, g, calcWt, antmap);

  xyz(1,:) = d.header.X(antmap);
  xyz(2,:) = d.header.Y(antmap);
  xyz(3,:) = d.header.Z(antmap);

  %------------------------------------------------------------
  % Finally, write the data!
  %------------------------------------------------------------

  climaxMatWriteUvf(uvffile, data, jd, srcname, coord, freq, xyz, double(uvw), uint32(length(antmap)), bw);

return

%------------------------------------------------------------
% Get the data in the form that climaxMatWriteUvf wants it
%------------------------------------------------------------

function [data uvw jds freq bw] = getData(d, nAnt, nBase, groups, calcWt, antmap)

  % Check for random inserts of groups in the miriad file.  

  if(mod(length(groups)/(nBase), 1) ~= 0)
    groups = purgeOutOfOrderGroups(groups, antmap);
  end

  nInt = length(groups)/(nBase);

  nFreq = size(groups(1).vis, 1);
  nChan = size(groups(1).vis, 2);
  flag  = [groups(:).flag];

  sfreq  = groups(1).sfreq;
  ischan = groups(1).ischan;
  sdf    = groups(1).sdf;

  centerfreq = sfreq;

  goodFreq = zeros(1,nFreq);

  iGood = 1;
  for iFreq=1:nFreq
    if(length(find(flag(iFreq,:)==0)) > 0)
      goodFreq(iFreq) = 1; 
      freqInds(iGood) = iFreq;
      iGood = iGood+1;
    end

    startchan = 1+(iFreq-1)*nChan
    stopchan  = startchan + nChan-1
    midchan   = (stopchan+startchan)/2
    centerfreq(iFreq) = sfreq(iFreq) + sdf(iFreq) * (midchan - double(ischan(iFreq)))

  end

  nGoodFreq = sum(goodFreq);
  vis = zeros(nInt, nBase, nGoodFreq, 3);

  vis = [groups(:).vis];
  vis = reshape(vis, nFreq, nChan, nBase, nInt);
  vis = permute(vis, [4 3 1 2]);

  flags = [groups(:).flag];
  flags = reshape(flags, nFreq, nChan, nBase, nInt);
  flags = permute(flags, [4 3 1 2]);

  vis(flags==0) = complex(nan,nan);
  visbar        = nanmean(vis,4);

  revar = nanvar(real(vis),[],4);
  imvar = nanvar(imag(vis),[],4);

  nchan = nansum(flags, 4);
  var = (revar + imvar)/2;
  varwt = var./nchan;

  temps = [groups(:).systemp];
  temps = reshape(temps, 2, nFreq, nBase, nInt);
  temps = permute(temps, [4 3 2 1]);

  ant1 = [groups(:).ant1];  
  ant1 = reshape(ant1, nBase, nInt);
  ant1 = permute(ant1, [2 1]);

  ant2 = [groups(:).ant2];  
  ant2 = reshape(ant2, nBase, nInt);
  ant2 = permute(ant2, [2 1]);

  jyperk = d.header.jyperka;

  jyfac2 = jyperk(ant1) .* jyperk(ant2);
  jyfac2 = repmat(jyfac2, [1,1,nFreq]);
  tsys2  = temps(:,:,:,1) .* temps(:,:,:,2);

  tint = [groups(:).inttime];  
  tint = reshape(tint, nBase, nInt);
  tint = permute(tint, [2 1]);
  tint = repmat(tint, [1,1,nFreq]);

  df = double(flags);
  df(flags==0) = nan;
  df = nansum(df,4);

  cf = repmat(sdf', nInt, nBase);
  cf = abs(reshape(cf, nInt, nBase, nFreq));
  
  % Calc bw in GHz

  bw = df .* cf * 1e9;
  bw = bw * 0.875;

  % Finally, put it all together

  calcvar = (tsys2 .* jyfac2) ./ (bw .* tint);
  calcvar = calcvar/2;
  
  if calcWt
    wt = 1.0/calcvar;
  else 
    wt = 1.0/varwt;
  end

  data(:,:,:,1) = real(visbar);
  data(:,:,:,2) = imag(visbar);
  data(:,:,:,3) = wt;

  u = [groups(:).u];
  u = reshape(u, nBase, nInt);
  u = permute(u, [2 1]);

  v = [groups(:).v];
  v = reshape(v, nBase, nInt);
  v = permute(v, [2 1]);

  uvw(:,:,1) = u*1e-9;
  uvw(:,:,2) = v*1e-9;
  uvw(:,:,3) = v*1e-9;

  mjds = [groups(:).mjd];
  mjds = reshape(mjds, nBase, nInt);
  mjds = permute(mjds, [2 1]);
  mjds = mjds(:,1);

  jd = mjds + 2400000.5;
  jds(:,1) = double(uint32(jd));
  jds(:,2) = double(jd - jds(:,1));

  data = double(data(:,:,freqInds,:));
  freq = centerfreq(freqInds)' * 1e9;

  df = double(flags);
  df(flags==0) = nan;
  df = nansum(df,4);
  df(df==0) = nan;
  df = nanmean(df,2);
  df = nanmean(df,1);

  meandf = nanmean(df);
  df(isnan(df)) = meandf;
  bw = squeeze(df) .* abs(squeeze(sdf*1e9));
  bw = bw';

  % Replace any zero bandwidth with the average of the non-zero values

  nonzeroinds = find(bw ~= 0);
  zeroinds = find(bw == 0);
  meanbw = mean(bw(nonzeroinds));
  bw(zeroinds) = meanbw;

return

%------------------------------------------------------------
% Return groups and coordinates corresponding to srcname
%------------------------------------------------------------

function [g coord] = getSrc(d, srcname)
  inds = [];

  nSrc = length(d.header.sources);

  for iSrc=1:nSrc
    if(strcmp(d.header.sources(iSrc).name, srcname))
      inds = d.header.sources(iSrc).inds;  					 
      coord(1) = d.header.sources(iSrc).obsra/pi * 12.0;
      coord(2) = d.header.sources(iSrc).obsdec/pi * 180.0;
      coord(3) = d.groups(1).mjd;
    end
  end

  if(length(inds) == 0)
    error(sprintf('No source %s found', srcname));
  end

  g = d.groups(inds);

  % Make sure this is sorted by time order, since UVF readers may 
  % expect it

  mjds = [g(:).mjd];
  [mjdssorted, I] = sort(mjds);
  g = [g(I)];

return

%------------------------------------------------------------
% Return various mappings for a subset of groups
%------------------------------------------------------------

function [nAnt antmap revmap bases] = getBaselineMapping(g)
  [nAnt antmap revmap] = getAntennaIndexArray(g);
%  groups = getGroupIndexArray(g);
  bases = getBaselineIndexArray(g);
return

%------------------------------------------------------------
% Get the antenna mappings
%------------------------------------------------------------

function [nAnt antmap revmap] = getAntennaIndexArray(g)

  ng = length(g);
  go = 1;
  iGroup = 1;
  nAnt = 0;

  ants1 = [g(:).ant1];
  ants2 = [g(:).ant2];
  ants = [ants1 ants2];
  nAnt = length(unique(ants));

  %------------------------------------------------------------
  % Unique returns an ordered array of unique CARMA-indexed
  % antenna numbers
  %------------------------------------------------------------

  antmap = unique(ants);

  for iAnt=1:nAnt
    revmap(antmap(iAnt)) = iAnt;
  end

return

function bases = getGroupIndexArray(g)
  [nAnt antmap reversemap] = getAntennaIndexArray(g);

  iBase=1;

  for iAnt1=1:nAnt
    ant1 = reversemap(iAnt1);
    ant2 = ant1;
    bases(ant1, ant2) = iBase;
    iBase = iBase+1;
  end

  for iAnt1=1:nAnt-1
    for iAnt2=iAnt1+1:nAnt
      ant1 = reversemap(iAnt1);
      ant2 = reversemap(iAnt2);
      bases(ant1, ant2) = iBase;
      iBase = iBase+1;
    end
  end

return

function bases = getBaselineIndexArray(g)
  [nAnt antmap reversemap] = getAntennaIndexArray(g);

  iBase=1;

  for iAnt1=1:nAnt-1
    for iAnt2=iAnt1+1:nAnt
      ant1 = antmap(iAnt1);
      ant2 = antmap(iAnt2);
      bases(ant1, ant2) = iBase;
      iBase = iBase+1;
    end
  end

return

function gret = purgeOutOfOrderGroups(g, antmap)

  nAnt = length(antmap);
  nBase = (nAnt*(nAnt-1))/2;

  iBase = 1;
  for iAnt1=1:nAnt-1
    for iAnt2=iAnt1+1:nAnt
      bases(iBase, 1) = antmap(iAnt1);
      bases(iBase, 2) = antmap(iAnt2);
      iBase = iBase+1;
    end
  end

  nGroup = length(g);
  nGood  = 0;

  for iGroup=1:nGroup

    iBase = mod(nGood, nBase)+1;
    group = g(iGroup);

    if(group.ant1 ~= bases(iBase,1) | group.ant2 ~= bases(iBase,2))
      sprintf('Ind %d is bad', iGroup)
    else
      nGood = nGood+1;
      goodInds(nGood) = iGroup;
    end

  end

  gret = g(goodInds);

return

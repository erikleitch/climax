function [data s1 s2] = readArcShadowTest(startUtc, stopUtc)

addpath matlab;

regs = {'array.frame.utc double',...
	'antenna*.tracker.actual double',...
	'array.carma.config',...
	'array.carma.pad',...
	'array.carma.nPad',...
	'array.carma.antType',...
        'array.sza.config',...
        'array.sza.pad'};

addpath /home/szadaq/sza_analysis/reduc

data = read_arc(startUtc, stopUtc, regs);

az = squeeze(data.antenna.tracker.actual(:,1,:));
el = squeeze(data.antenna.tracker.actual(:,2,:));

szaMatCarmaConfiguration(data);

keyboard

s1 = szaMatShadowCarmaNew(az, el, data);
s2 = szaMatShadowCarma(az, el, 'I', 'DO');

return




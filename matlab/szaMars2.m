function [flux error]=szaMars2(type, mjd, freq, marsDir, marsFile)

if(~exist('marsDir'))
  marsDir=[];
end

if(~exist('marsFile'))
  marsFile=[];
end

if(isempty(marsDir))
  marsDir='/home/szadaq/sza/array/ephem/';
end

if(isempty(marsFile))
  marsFile='mars_new.mod';
end

% Call mex function

[flux error] = szaMatMars(marsDir, marsFile, type, mjd, freq);

return

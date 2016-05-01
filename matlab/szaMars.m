function [flux error]=szaMars(type, mjd, freq, marsDir, marsFile)

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
%  marsFile='mars.mod';
  marsFile='mars_new_bb.mod';
end

% Call mex function

%[flux error] = szaMatMars(marsDir, marsFile, type, mjd, freq);
[flux error] = szaMatMarsNew(marsDir, marsFile, type, mjd, freq);

return

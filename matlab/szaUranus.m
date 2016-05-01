function [flux error]=szaUranus(type, mjd, freq, urnsDir, urnsFile)

if(~exist('urnsDir'))
  urnsDir=[];
end

if(~exist('urnsFile'))
  urnsFile=[];
end

if(isempty(urnsDir))
  urnsDir='/home/szadaq/sza/array/ephem/';
end

if(isempty(urnsFile))
  urnsFile='uranus_ediam.mod';
end

% Call mex function

[flux error] = szaMatUranus(urnsDir, urnsFile, type, mjd, freq);

return

function data=szaReadArc(regs,start,finish,arcdir,calfile)
% data=szaReadArc(regs,start,end,arcdir,calfile)
%
% Extract a given set of registers between given start and stop times
% from arc files into matlab structure
%
% regs=cell array of register names of form board.reg[x-y]
% single indecies or range allowed. Note that these are archive native
% zero based indecies, whereas returned matlab data has 1 based indecies.
% start/finish=start finish times as strings with format e.g.:
% 01-Jan-2005:00:00:00
% optional arcdir=directory containing archive files
% optional calfile=cal file
%
% e.g.: regs={'array.frame.utc','array.weather.airTemperature'}
%       d=szaReadArc(regs,'25-Dec-2004:00:00:00','26-Dec-2004:07:00:00')
%       plot(d.array.frame.utc,d.array.weather.airTemperature)

if(~exist('arcdir'))
  arcdir=[];
end
if(~exist('calfile'))
  calfile=[];
end
if(isempty(arcdir))
  arcdir='/data/szadaq/arc';
end
if(isempty(calfile))
  calfile='/home/szadaq/sza/array/conf/cal';
end

% Ensure regs unique
if(any(size(unique(regs))~=size(regs)))
  error('regs should be unique');
end

% Call mex function
data=szaMatReadArc(regs,start,finish,arcdir,calfile);

return

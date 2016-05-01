function [ra dec] = app2j2000(radeg, decdeg, mjd)
  [ra dec] = climaxMatAppToJ2000(radeg/180*12, decdeg, mjd);
  ra = ra/12*180;
return

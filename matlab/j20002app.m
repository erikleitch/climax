function [ra dec] = j20002app(radeg, decdeg, mjd)
  [ra dec] = climaxMatJ2000ToApp(radeg/180*12, decdeg, mjd);
  ra = ra/12*180;
return

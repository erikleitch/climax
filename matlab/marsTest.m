mt = dlmread('/home/szadaq/sza/array/ephem/mars_test.mod');

plot(mt(:,1), mt(:,4))

[t26, e26] = szaMars3('temp', mt(:,1), 26);
[t32, e32] = szaMars3('temp', mt(:,1), 32);

hold on

plot(mt(:,1), t26, 'r')
plot(mt(:,1), t32, 'b')

%plot(mt(:,1), mt(:,4)-t26, 'c')

mjdint = min(mt(:,1)):0.1:max(mt(:,1));

[t26int, e26int] = szaMars3('temp', mt(:,1), 26);
[t29int, e29int] = szaMars3('temp', mt(:,1), 29);

[t86int, e86int] = szaMars3('temp', mt(:,1), 86);

plot(mt(:,1), t29int, 'g')
plot(mt(:,1), t86int, 'm')

[t100int, e100int] = szaMars3('temp', mt(:,1), 100);

plot(mt(:,1), t100int, 'y')

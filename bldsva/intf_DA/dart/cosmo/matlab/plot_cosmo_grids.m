fname = '../data/lfff00030000.nc';


rlon = ncread(fname,'rlon');
rlat = ncread(fname,'rlat');
[rx, ry] = meshgrid(rlon,rlat);

subplot(2,1,1)
plot(rx(:), ry(:), 'x');
title('rotated')
xlabel('rlon')
ylabel('rlat')

lon = ncread(fname,'lon');
lat = ncread(fname,'lat');

subplot(2,1,2)
plot(lon(:), lat(:), 'o');
title('geographic')
xlabel('lon')
ylabel('lat')


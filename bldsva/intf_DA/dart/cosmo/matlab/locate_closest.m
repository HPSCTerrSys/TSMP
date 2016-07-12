cosmogrid = load('../work/exhaustive_grid.m');

indI = cosmogrid(:,1);
indJ = cosmogrid(:,2);
indK = cosmogrid(:,3);

indiii = cosmogrid(:,4);
ulon   = cosmogrid(:,5);
ulat   = cosmogrid(:,6);
ulev   = cosmogrid(:,7);
rlon   = cosmogrid(:,8);
rlat   = cosmogrid(:,9);

% pertains to k = 2
index1 = 601;
indexN = 1200;

figure(1); orient landscape

plot(ulon(index1:indexN),ulat(index1:indexN),'x')
title('geographic grid for U stagger')
xlabel('longitude')
ylabel('latitude')

for i = 601:1200
   h = text(ulon(i),ulat(i),sprintf('%d',i));
   set(h,'HorizontalAlignment','center','VerticalAlignment','bottom')
end

grid

figure(2); orient landscape

plot(rlon(index1:indexN),rlat(index1:indexN),'x')
title('rotated grid for U stagger')
xlabel('rotated longitude')
ylabel('rotated latitude')

for i = 601:1200
   h = text(rlon(i),rlat(i),sprintf('%d',i));
   set(h,'HorizontalAlignment','center','VerticalAlignment','bottom')
end

grid


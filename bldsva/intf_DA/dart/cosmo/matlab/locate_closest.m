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

plot(ulon(index1:indexN),ulat(index1:indexN),'o')

for i = 601:1200
   text(ulon(i),ulat(i),sprintf('%d',i));
end



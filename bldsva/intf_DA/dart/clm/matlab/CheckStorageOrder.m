function bob = CheckStorageOrder(fname,varname)
%% DART CheckStorageOrder - explore xxx(column,level) variable order.
%
% fname = 'clm_restart.nc';
% varname = 'H2OSOI_LIQ';
% x = CheckStorageOrder(fname,varname);

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id: CheckStorageOrder.m 6265 2013-06-13 19:38:46Z thoar $

sparselat = nc_varget(fname,'cols1d_lat');
sparselon = nc_varget(fname,'cols1d_lon');
datmat    = nc_varget(fname,varname);
datmat(datmat > 1.0E30) = NaN;

nlev = size(datmat,2);

for i = 1:6
   figure(i);
   plot3(sparselon,sparselat,datmat(:,i),'x');
   h = title(sprintf('%s fortran/matlab index %d of %d',varname,i,nlev));
   set(h,'Interpreter','none');
   hold on; worldmap; grid; view(0,90)
end

% <next few lines under version control, do not edit>
% $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/clm/matlab/CheckStorageOrder.m $
% $Revision: 6265 $
% $Date: 2013-06-13 21:38:46 +0200 (Thu, 13 Jun 2013) $


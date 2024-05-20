function narray = myArrayResample3D(oldpixsize,array,newpixsize,intmethod)
% nimg = imresample(oldpixsize,img,newpixsize,intmethod)
% This function resamples the images at the new grid points
% defined by the new pixel sizes. It assumes that intensities are
% defined at the pixel centers
%
%    array  : original 3D array to be resampled
%    narray : newly sampled image
%    oldpixsize : a vector of the form [xpixsize ypixsize zpixsize] 
%       for the original image, e.g., [0.5,0.5 0.5]
%    newpixsize : is a vector of the form [xpixsize ypixsize zpixsize]
%       for the new image, e.g., [0.2,0.2 0.2]
%    intmethod: same as interp2
%       'nearest' - nearest neighbor
%       'linear'  - bilinear
%       'cubic'   - bicubic
%       'spline'  - spline

% Omer Demirkaya 12/14/2008
% modified by Susanne Rafelski 2/20/2009

% Find the dimensions of the image
[r,c,d] = size(array);
%r = r-1; c = c-1; d = d-1;

% smaller variable names
ops = oldpixsize;
nps = newpixsize;
ss  = ops/2;
nn  = nps/2;
% create the meshes for old and new image

%let meshes in x,y (ie where size is reduced) go from half-pix to size-halfpix
%let mesh in z (ie where size is increased) go from 1 to size
[Ox,Oy, Oz] = meshgrid(ss(1):ops(1):c*ops(1)-ss(1),...
    ss(2):ops(2):r*ops(2)-ss(2), ...
    1:ops(3):d);
[Nx,Ny, Nz] = meshgrid(nn(1):nps(1):c*ops(1)-nn(1),...
    nn(2):nps(2):r*ops(2)-nn(2), ...
    1:nps(3):d);
% create the new image
narray = interp3(Ox,Oy, Oz ,array,Nx,Ny, Nz,intmethod);
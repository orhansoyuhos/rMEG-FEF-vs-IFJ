function tutorial_nwa_connectivityviewer_save_figures(data, param, colorlim, seed_target, flag, pos2d)

if nargin<2,
  param = 'cohspctrm';
end
if nargin<3,
  tmp = data.(param);
  if all(tmp(:)>=0), 
    colorlim = [0 max(tmp(:))];
  else
    colorlim = [-1 1].*max(abs(tmp(:)));
  end
end

if isfield(data, 'brainordinate') && ~isfield(data, 'pos')
  % this seems to be parcellated data, 'unparcellate' on the fly
  tmpdat = data.(param);
  data   = data.brainordinate;
  atlas  = data;
  
  dat = zeros(size(data.pos,1));
  for k = 1:numel(data.parcellationlabel)
    for m = 1:numel(data.parcellationlabel)
      dat(data.parcellation==k,data.parcellation==m) = tmpdat(k,m);
    end
  end
else
  dat = data.(param);
  atlas = [];
end

% create sphere for location highlighting
[x,y,z] = sphere(10);

figure; hold on;
ft_plot_mesh(data,'vertexcolor', mean(dat,2))
h1 = gca;
h2 = get(h1, 'children');
if flag == false
    h2.EdgeColor = [0.65,0.65,0.65]; %[0.80,0.80,0.80];%[0.94,0.94,0.94];
    h2.LineStyle = ':';
end
h3 = surf(h1, 2*x,2*y,2*z+100, 'facecolor', 'w', 'edgecolor', 'none');

% % Orhan:
% change the view angle
if seed_target(1) == 'r'
    view(360,0)
elseif seed_target(1) == 'l'
    view(180,0)
end 

% Orhan:
cbar = colorbar;
cbar.FontSize = 14;
set(cbar.XLabel,{'String','Rotation','Position','FontWeight'},{'Z',0,[0.5 colorlim(1)-0.1],'bold'})

% % Orhan:
fancyplotting([h1 h2 h3], data.pos, data.tri, dat, atlas, seed_target, flag, pos2d);

caxis(colorlim);
end

function fancyplotting(axh, pos, tri, dat, atlas, seed_target, flag, pos2d)

% get the clicked position in axis coordinates
% pos2d = get(axh(1),'CurrentPoint');

% get the intersection with the mesh
[ipos, d] = intersect_line(pos, tri, pos2d(1,:), pos2d(2,:));
[md, ix]  = min(abs(d));

% get the closest mesh point, and the corresponding index
dpos = pos - ipos(ix*ones(size(pos,1),1),:);
indx = nearest(sum(dpos.^2,2),0);

x = get(axh(3),'xdata');
x = x-mean(x(:));
y = get(axh(3),'ydata');
y = y-mean(y(:));
z = get(axh(3),'zdata');
z = z-mean(z(:));

set(axh(3),'xdata', x+pos(indx,1));
set(axh(3),'ydata', y+pos(indx,2));
set(axh(3),'zdata', z+pos(indx,3));

% update the functional data
set(axh(2), 'facevertexcdata', dat(:, indx));

% Orhan: Not used anymore, but keep it for later
% for i = 1:numel(atlas.parcellationlabel)
%     scout_idx = (atlas.parcellation == i);
%     P = atlas.pos(scout_idx,:);
%     k = boundary(P);
%     hold on;
%     trisurf(k,P(:,1),P(:,2),P(:,3),'Facecolor','black','FaceAlpha', 0.1,'EdgeColor', 'none')
% end

% Orhan: Only for ROIs' contours; you need add 'ROI_idx' parameter to func
% for i = 1:length(ROI_idx)
%     scout_idx = (atlas.parcellation == ROI_idx(i));
%     ft_plot_mesh(atlas, 'vertexcolor', dat(:, indx), 'contour', scout_idx, 'contourcolor', 'k', 'contourlinewidth', 0.5);
% end

% % Orhan: for contour of the parcels
if flag == true
    if or(strcmp(seed_target, 'right-right'), strcmp(seed_target, 'left-right'))
        for i = 1:180
            mask = atlas.parcellation == i+180;
            scout_idx{i,:} = double(mask);
        end
        ft_plot_mesh(atlas, 'vertexcolor', dat(:, indx), 'contour', scout_idx, 'contourcolor', [185, 185, 185]/255, 'contourlinewidth', 0.5);
    else
        for i = 1:180
            mask = atlas.parcellation == i;
            scout_idx{i,:} = double(mask);
        end
        ft_plot_mesh(atlas, 'vertexcolor', dat(:, indx), 'contour', scout_idx, 'contourcolor', [185, 185, 185]/255, 'contourlinewidth', 0.5);
    end
end

if isempty(atlas)
  fprintf('you clicked on position %d with coordinates %3.2f %3.2f %3.2f\n', [indx pos(indx,:)]);
else
  fprintf('you clicked on parcel %s\n', atlas.parcellationlabel{atlas.parcellation(indx)});
  title(atlas.parcellationlabel{atlas.parcellation(indx)}, 'Interpreter', 'none');
end





end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the code below is copied from fieldtrip/private

function [points, pos, indx] = intersect_line(pnt, tri, pnt1, pnt2)

% INTERSECT_LINE finds the intersection points between a mesh and a line.
%
% Use as:
%   [points, pos, indx] = intersect_line(pnt, tri, pnt1, pnt2)
% 
% Where pnt (Nx3) and tri (Mx3) define the mesh, and pnt1 (1x3) and pnt2
% (1x3) define the line. The output argument points (Px3) are the
% intersection points, pos (Px1) the location on the line (relative to
% pnt1) and indx is the index to the triangles of the mesh that are
% intersected.
%
% This code is based from a function from the geom3d toolbox, that can be
% found on matlab's file exchange. The original help is pasted below. The
% original function was released under the BSD-license.
%
% Adapted to FieldTrip by Jan-Mathijs Schoffelen 2012

%INTERSECTLINEMESH3D Intersection points of a 3D line with a mesh
%
%   INTERS = intersectLineMesh3d(LINE, VERTICES, FACES)
%   Compute the intersection points between a 3D line and a 3D mesh defined
%   by vertices and faces.
%
%   [INTERS POS INDS] = intersectLineMesh3d(LINE, VERTICES, FACES)
%   Also returns the position of each intersection point on the input line,
%   and the index of the intersected faces.
%   If POS > 0, the point is also on the ray corresponding to the line. 
%   
%   Example
%   intersectLineMesh3d
%
%   See also
%   meshes3d, triangulateFaces
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-12-20,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


tol = 1e-12;

ntri = size(tri,1);

% normals to the triangles
n   = normals(pnt, tri, 'triangle');

% vectors describing the edges
t0  = pnt(tri(:,1),:);
u   = pnt(tri(:,2),:) - t0;
v   = pnt(tri(:,3),:) - t0;

% get the direction of the line
d   = pnt2 - pnt1;
d   = d/norm(d);

% vector between triangle origin and line origin
w0  = pnt1(ones(ntri,1),:) - t0;

a   = -sum(n.*w0, 2);
b   =  n*d';

valid = abs(b)>tol & sqrt(sum(n.^2,2))>tol;

% compute intersection point of line with supporting plane
% If pos < 0: point before ray
% IF pos > |dir|: point after edge
pos = a ./ b;

% coordinates of intersection point
points = pnt1(ones(ntri,1),:) + pos*d;

% normalize direction vectors of triangle edges
uu  = sum(u.^2, 2);
uv  = sum(u.*v, 2);
vv  = sum(v.^2, 2);

% coordinates of vector v in triangle basis
w   = points - t0;
wu  = sum(w.*u, 2);
wv  = sum(w.*v, 2);

% normalization constant
D = uv.^2 - uu .* vv;

% test first coordinate
s    = (uv .* wv - vv .* wu) ./ D;
ind1 = s < 0.0 | s > 1.0;
points(ind1, :) = NaN;
pos(ind1)       = NaN;

% test second coordinate, and third triangle edge
t    = (uv .* wu - uu .* wv) ./ D;
ind2 = t < 0.0 | (s + t) > 1.0;
points(ind2, :) = NaN;
pos(ind2)       = NaN;

% keep only interesting points
inds   = ~ind1 & ~ind2 & valid;
points = points(inds, :);

pos  = pos(inds);
indx = find(inds);

end

function [nrm] = normals(pnt, tri, opt)

% NORMALS compute the surface normals of a triangular mesh
% for each triangle or for each vertex
%
% [nrm] = normals(pnt, tri, opt)
% where opt is either 'vertex' or 'triangle'

% Copyright (C) 2002-2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<3
  opt='vertex';
elseif (opt(1)=='v' | opt(1)=='V')
  opt='vertex';
elseif (opt(1)=='t' | opt(1)=='T')
  opt='triangle';
else
  error('invalid optional argument');
end

npnt = size(pnt,1);
ntri = size(tri,1);

% shift to center
pnt(:,1) = pnt(:,1)-mean(pnt(:,1),1);
pnt(:,2) = pnt(:,2)-mean(pnt(:,2),1);
pnt(:,3) = pnt(:,3)-mean(pnt(:,3),1);

% compute triangle normals
% nrm_tri = zeros(ntri, 3);
% for i=1:ntri
%   v2 = pnt(tri(i,2),:) - pnt(tri(i,1),:);
%   v3 = pnt(tri(i,3),:) - pnt(tri(i,1),:);
%   nrm_tri(i,:) = cross(v2, v3);
% end

% vectorized version of the previous part
v2 = pnt(tri(:,2),:) - pnt(tri(:,1),:);
v3 = pnt(tri(:,3),:) - pnt(tri(:,1),:);
nrm_tri = cross(v2, v3);


if strcmp(opt, 'vertex')
  % compute vertex normals
  nrm_pnt = zeros(npnt, 3);
  for i=1:ntri
    nrm_pnt(tri(i,1),:) = nrm_pnt(tri(i,1),:) + nrm_tri(i,:);
    nrm_pnt(tri(i,2),:) = nrm_pnt(tri(i,2),:) + nrm_tri(i,:);
    nrm_pnt(tri(i,3),:) = nrm_pnt(tri(i,3),:) + nrm_tri(i,:);
  end
  % normalise the direction vectors to have length one
  nrm = nrm_pnt ./ (sqrt(sum(nrm_pnt.^2, 2)) * ones(1,3));
else
  % normalise the direction vectors to have length one
  nrm = nrm_tri ./ (sqrt(sum(nrm_tri.^2, 2)) * ones(1,3));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product to replace the MATLAB standard version
function [c] = cross(a,b)
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end

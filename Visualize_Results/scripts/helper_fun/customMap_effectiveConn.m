function CMap = customMap_effectiveConn(N)

% % Orhan:
% https://jdherman.github.io/colormap/

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Sylvain Baillet

if (nargin < 1)
    N = 128; % Number of color levels
end

%%

Wwidth = 0;  % Width of the white zero level
Half = floor(N/2);

R1  = linspace(0.1569, 1, Half-Wwidth);
R2  = ones(1,Wwidth*2);
R3  = flip(linspace(0, 1, Half-Wwidth));
R   = [R1 R2 R3]';

G1  = linspace(0, 1, Half-Wwidth);
G2  = ones(1,Wwidth*2);
G3  = flip(linspace(0.3137, 1, Half-Wwidth));
G   = [G1 G2 G3]';

B1  = linspace(0.3137, 1, Half-Wwidth);
B2  = ones(1,Wwidth*2);
B3  = flip(linspace(0.1569, 1, Half-Wwidth));
B   = [B1 B2 B3]';

CMap = [R G B];



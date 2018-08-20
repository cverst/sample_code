function [faces, verts] = basePlate(S,thr,rsfact,scfact,rotation,position,flag,FN)
% BASEPLATE - constructs base plate for in vitro cochlear preparation
%   [faces, vertices] = basePlate(S,thr,rsfact,scfact,rotation,position,flag,FN)
%   constructs a base plate to be printed by a 3D printer. The
%   negative shape of a structure is cut out by thresholding a series of
%   images preprocessed by dicomview.
%
%   INPUTS:   S        - output struct from dicomview
%             thr      - threshold used as edge of structure
%                        [default: threshold from dicomview output]
%             rsfact   - resample factor of original image stack resolution
%                        [default: 1]
%             scfact   - scale factor to adjust size of cutout [default: 1]
%             rotation - vector specifying x, y, and z rotation [x y z]
%                        [default: [0 0 0]]
%             position - translation across baseplate of rotated object
%                        [default: [0 0 0]]
%             flag     - can be (a combination of) '-flip', to invert the
%                        top and bottom of the rotated object, and '-small',
%                        to reduce the dimensions of the baseplate for
%                        faster evaluation of the function
%                        [default: '-none']
%             FN       - filename of output .stl file
%                        [default: 'test']
%
%   OUTPUTS: none     - writes an .stl file of the base plate and a second
%                       .stl file for the inverted shape (other ear) to the
%                       current directory
%            faces    - faces of the base plate shape
%            vertices - vertices of the base plate shape
%
%   See also dicomview, isosurface, stlwrite.
%

% Defaults
if nargin < 2 || isempty(thr), thr = S(1).contourVal; end
if nargin < 3 || isempty(rsfact), rsfact = 1; end
if nargin < 4 || isempty(scfact), scfact = 1; end
if nargin < 5 || isempty(rotation), rotation = [0 0 0]; end
if nargin < 6 || isempty(position), position = [0 0 0]; end
if nargin < 7 || isempty(flag), flag = '-none'; end
if nargin < 8 || isempty(FN), FN = 'test'; end

% Timer
tic

% Construct 3D array
S1 = S(1);
rows = double(S1.info.Rows);
cols = double(S1.info.Columns);
stacks = numel(S);
MM = reshape([S.imageAdj],rows,cols,stacks);
MM = double(MM);
MM(MM<thr) = 0; % threshold matrix to binary
MM(MM>=thr) = 1; % threshold matrix to binary

% Scale factors
Xscale = S1.info.PixelSpacing(1)*1e3*scfact;
Yscale = S1.info.PixelSpacing(2)*1e3*scfact;
Zscale = S1.info.SliceThickness*1e3*scfact;
clear S % memory issue

% Timer
disp('Build 3D array')
toc

% Rotate array
if any(rotation~=0)
    MM = local_rotate(MM,rotation);
end
if any(strfind(flag,'-flip'))
    MM = MM(:,:,numel(MM(1,1,:)):-1:1);
end

% Timer
disp('Rotation')
toc

% Make axes
sz = size(MM);
Xax = (0:sz(2)-1)*Yscale;
Yax = (0:sz(1)-1)*Xscale;
Zax = (0:sz(3)-1)*Zscale;

% Create grid and resample matrix if queried
[MM, Xax, Yax, Zax] = local_gridresample(MM,Xax,Yax,Zax,rsfact);

% Timer
disp('Grid/resampling')
toc

% Create base plate with cutout
[BP, XX, YY, ZZ] = local_base(MM,Xax,Yax,Zax,position.*[1 1 scfact],flag);

% Timer
disp('Base plate')
toc

% Compute surfaces
FV = isosurface(XX,YY,ZZ,BP,.9);
FVmirror = FV;
FVmirror.vertices(:,:) = FVmirror.vertices.*repmat([-1,1,1],numel(FVmirror.vertices(:,1)),1);

% Timer
disp('Isosurface')
toc

% Output if necessary
if nargout > 0
    faces = FV.faces;
    verts = FV.vertices;
else
    stlwrite([FN '.stl'],FV.faces,FV.vertices)
    stlwrite([FN '_mirror.stl'],FVmirror.faces,FVmirror.vertices)
    % Timer
    disp('Write file')
    toc
end


function [MM, Xax, Yax, Zax] = local_gridresample(MM,Xax,Yax,Zax,rsfact)
% Create meshgrid and resample data if necessary

% Resample matrix if factor ~= 1
if rsfact ~= 1
    [XX, YY, ZZ] = meshgrid(Xax,Yax,Zax);
    nrows = numel(Yax);
    ncols = numel(Xax);
    nstacks = numel(Zax);
    Xax = linspace(min(Xax),max(Xax),round(ncols/rsfact));
    Yax = linspace(min(Yax),max(Yax),round(nrows/rsfact));
    Zax = linspace(min(Zax),max(Zax),round(nstacks/rsfact));
    [XXnew, YYnew, ZZnew] = meshgrid(Xax,Yax,Zax);
    MM = interp3(XX,YY,ZZ,MM,XXnew,YYnew,ZZnew);
end


function [BP, XX, YY, ZZ] = local_base(MM,Xax,Yax,Zax,pos,flag)
% Creates plate and with cochlea cutout in matrix form

% Plate dimensions
PH = .8; % height (mm)
if any(strfind(flag,'-small'))
    PL = 1; % length (mm)
    PW = 1; % width (mm)
else
    PL = 18; % length (mm)
    PW = 10; % width (mm)
end

% Grid spacing
dx = diff(Xax(1:2)); % um
dy = diff(Yax(1:2)); % um
dz = diff(Zax(1:2)); % um

% Plate dimensions in grid indices
[nrows, ncols, junk] = size(MM); %#ok<NASGU>
ibpheight = round(pos(3)*1e3/dz)+1; % starting index of base plate
nLind = max([ncols round(PL*1e3/dx)]); % number of indices LENGTH / X
nWind = max([nrows round(PW*1e3/dy)]); % number of indices WIDTH / Y
nHind = round(PH*1e3/dz); % number of indices HEIGHT / Z

% Create array containing plate
BP = zeros(nWind+2,nLind+2,nHind+2); % create new matrix
BP(2:end-1,2:end-1,2:end-1) = 1;

% Put cochlea into base plate matrix
rowindx = (1:nrows)+round(pos(2)*1e3/dy);
colindx = (1:ncols)+round(pos(1)*1e3/dx);
rowindx(rowindx>nWind+2) = []; % cut off values outside base plate
colindx(colindx>nLind+2) = []; % cut off values outside base plate
BP(rowindx,colindx,2:end-1) = BP(rowindx,colindx,2:end-1) + MM(1:numel(rowindx),1:numel(colindx),ibpheight:ibpheight+nHind-1)*2;

% New meshgrid
[XX, YY, ZZ] = meshgrid((0:nLind+1)*dx,(0:nWind+1)*dy,(0:nHind+1)*dz);

% Keep base plate only
BP(BP~=1) = 0;

% Find largest structure and remove others
CC = bwconncomp(BP,18); % compute connectivity
numPixels = cellfun(@numel,CC.PixelIdxList);
[junk, indx] = max(numPixels); %#ok<ASGLU>
BP = zeros(size(BP));
BP(CC.PixelIdxList{indx}) = 1;


function Mfinal = local_rotate(M,rot)
%
% http://blogs.mathworks.com/steve/2006/08/17/spatial-transformations-three-dimensional-rotation/
%

% Make tfrom
% Rotation
xr = rot(1);
yr = rot(2);
zr = rot(3);
Tx = [1         0         0         0
      0         cosd(xr)  sind(xr)  0
      0         -sind(xr) cosd(xr)  0
      0         0         0         1];
Ty = [cosd(yr)  0         -sind(yr) 0
      0         1         0         0
      sind(yr)  0         cosd(yr)  0
      0         0         0         1];
Tz = [cosd(zr)  sind(zr)  0         0
      -sind(zr) cosd(zr)  0         0
      0         0         1         0
      0         0         0         1];
% Translation
[rows, cols, stacks] = size(M);
[XX, YY, ZZ] = meshgrid(1:cols,1:rows,1:stacks);
cX = mean(XX(M==1));
cY = mean(YY(M==1));
cZ = mean(ZZ(M==1));
T1 = [1   0   0   0
      0   1   0   0
      0   0   1   0
      -cX -cY -cZ 1];
T2 = [1   0   0   0
      0   1   0   0
      0   0   1   0
      cX  cY  cZ  1];
%
T = T1*Tx*Ty*Tz*T2;
tform = maketform('affine',T);

% Transform
R = makeresampler('nearest', 'fill');
% Matrix bounds {origin xend yend zend xy xz yz xyz}
[nrows, ncols, nstacks] = size(M);
bnds = {[0;0;0;1] [ncols-1;0;0;1] [0;nrows-1;0;1] [0;0;nstacks-1;1] [ncols-1;nrows-1;0;1] [ncols-1;0;nstacks-1;1] [0;nrows-1;nstacks-1;1] [ncols-1;nrows-1;nstacks-1;1]};
newbnds = cellfun(@(V) ceil(T*V),bnds,'UniformOutput',false);
newbnds = [newbnds{:}].';
tsize_b = [max(newbnds(:,2))-min(newbnds(:,2)) max(newbnds(:,1))-min(newbnds(:,1)) max(newbnds(:,3))-min(newbnds(:,3))]+1;
%
F = 0;
Mrotated = tformarray(M,tform,R,[1 2 3],[1 2 3],tsize_b,[],F);

% Remove any empty planes
ikeepAll = any(Mrotated,3);
ikeepCol = any(ikeepAll,1);
ikeepRow = any(ikeepAll,2);
ikeepStack = any(squeeze(any(Mrotated,1)),1);
Mfinal = Mrotated(ikeepRow,ikeepCol,ikeepStack);

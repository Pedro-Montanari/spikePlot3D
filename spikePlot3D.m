function [Xdomain,Ydomain,V] = spikePlot3D(S,r0,gridS,gridD,varargin)
%spikePlot3D Generate spike-like 3D graphs in a structured grid from sample
% points in an unstructured grid
% 
%   This function takes sample points (defined by "x" and "y" coordinates
%   and a value "v") on an unstructured cartesian domain, distributes the 
%   value "v" of each sample point to a polar domain of radius "r0" 
%   according to a decaying exponential function, and sums the 
%   contributions of all polar domains to a structured cartesian domain.
%
%   The exponential function used in this algorithm is given by 
%   $$V = v*\left(1- \frac{r}{r_0}\right)\left(\frac{\exp(r_0-r)}{\exp(r_0)}\right)$$,
%   where "V" is the calculated value in the current point, "v" is the 
%   function value at the sample point, and "r" is the distance from the 
%   sample point to the current point. At the sample point the value of the
%   function is "v", and at the radius "r0" the value of the function is
%   zero.
%
%   Note: The purpose of this function is data visualization, and not
%         necessarily accuracy. The polar domains are approximately of 
%         radius "r0".
%
%   Usage:
%       >>  [Xdomain,Ydomain,V] = spikePlot3D(S,r0,gridS,gridD); % without
%       domain size explicitly defined
%       >>  [Xdomain,Ydomain,V] = spikePlot3D(S,r0,gridS,gridD,D); % with
%       domain size explicitly defined
%
%   Input arguments:
%       S : n-by-3 array of sample points. The 1st column is the
%           x-coordinate of the samples, the 2nd is the y-coordinate,
%           and the 3rd is the value. "n" is the number of sample points.
%       r0 : radius of the polar domains defined from the sample points.      
%       gridS : discretization of the polar domain. It is composed of
%               gridS-by-gridS points.
%       gridD : discretization of the structured cartesian domain. 
%               The function values of the polar domains are interpolated to
%               gridD-by-gridD points. For "n" sample points, 
%               (n*gridS)-by-(n*gridS) are interpolated to gridD-by-gridD 
%               points. To avoid losing resolution, set gridD at least 
%               equal to n*gridS.
%       D (OPTIONAL): Structured domain bounds, defined as a 2-by-2 array 
%                     in the form [Xmin Ymin; Xmax Ymax]. If no domain 
%                     bounds are specified, the resulting domain is the 
%                     bounding box that encompassess all the polar domains.
%   Output arguments:
%       Xdomain : gridD-by-gridD meshgrid array of the X coordinates
%                 of the structured domain
%       Ydomain : gridD-by-gridD meshgrid array of the Y coordinates
%                 of the structured domain
%       V : gridD-by-gridD meshgrid array of the V function values
%                 on the structured domain

% Unpack variables
Xos = S(:,1);
Yos = S(:,2);
Vs = S(:,3);
switch nargin
    case 5
        D = varargin{end}; % Domain limits
    case 4
    otherwise
        error("Wrong number of input arguments. Please specify 4 or 5");
end

% Create a unit circle with "gridS" points centered at the origin
[Xc,Yc,~] = sphere(gridS);
npts = height(S);

% For each sample point, create a circle with origin in the sample point.
% The "X" and "Y" coordinates of the polar domain of each sample are 
% obtained by scaling the unit circle to the radius "r0" and translating it
% to the origin of each sample point "Xos", "Yos". The Z coordinate is 
% calculated by the exponential function.
X = cell(npts,1);
Y = cell(npts,1);
Z = cell(npts,1);
x = cell(npts,1);
y = cell(npts,1);
z = cell(npts,1);

noDomainSpecified = false;
if ~exist("D","var") % If no domain bounds were specified
    noDomainSpecified = true;
    xmin = zeros(npts,1);
    xmax = zeros(npts,1);
    ymin = zeros(npts,1);
    ymax = zeros(npts,1);
end

for ii = 1:npts
    X{ii} = Xc*r0 + Xos(ii);
    Y{ii} = Yc*r0 + Yos(ii);
    r = sqrt((X{ii}-Xos(ii)).^2+(Y{ii}-Yos(ii)).^2);
    Z{ii} = Vs(ii)*(1-r./r0).*(exp(r0-r))/exp(r0);

    % Define the polar domain points as arrays. The domain points are
    % defined in a cartesian coordinate system
    x{ii} = reshape(X{ii},1,[])';
    y{ii} = reshape(Y{ii},1,[])';
    z{ii} = reshape(Z{ii},1,[])';

    % If domain coordinates were not specified, calculate them based on the
    % bounding box encompassing all polar domains.
    if noDomainSpecified
        xmin(ii) = min(X{ii},[],"all");
        xmax(ii) = max(X{ii},[],"all");
        ymin(ii) = min(Y{ii},[],"all");
        ymax(ii) = max(Y{ii},[],"all");
    end
end

% Create the structured domain
if noDomainSpecified % create the domain as the bounding box encompassing 
    % all polar domains.
    xdomain = linspace(min(xmin),max(xmax),gridD);
    ydomain = linspace(min(ymin),max(ymax),gridD);
else % use the bounds provided
    xdomain = linspace(D(1,1),D(2,1),gridD);
    ydomain = linspace(D(1,2),D(2,2),gridD);
end

% Initialize an empty domain grid and "n" temporary domain grids. The
% temporary domain grids contain the contribution of each sample point
% individually to the domain. The grids are added to the main domain grid
% to account for the contribution of all sample points together.
[Xdomain,Ydomain] = meshgrid(xdomain,ydomain);
V = zeros(length(Xdomain),height(Xdomain));
V_temp = cell(npts,1);

for ii = 1:npts
    % Interpolate the contribution of each individual sample to its
    % temporary domain grid.
    F = scatteredInterpolant(x{ii},y{ii},z{ii});
    V_temp{ii} = F(Xdomain,Ydomain);

    % Force the lower bound to zero
    V_temp{ii}(V_temp{ii}<0) = 0;

    % For each temporary grid, map the values to vary linearly
    % between zero and the value at the sample point
    V_temp{ii} = reshape(V_temp{ii},[1 numel(V_temp{ii})]);
    if noDomainSpecified
        V_temp{ii} = linmap(V_temp{ii},[0, max(Z{ii},[],"all")]);
    else
        V_temp{ii} = linmap(V_temp{ii},[0, max(V_temp{ii})]);
    end

    V_temp{ii} = reshape(V_temp{ii},[length(Xdomain),height(Xdomain)]);

    % Add the contribution of each temporary grid to the main domain grid
    V = V + V_temp{ii};
end

    function vout = linmap(vin,rout)
        % Source : https://machinelearning1.wordpress.com/2014/07/13/linear-vector-mapping-scaling-matlab/
        % function for linear mapping between two ranges
        % inputs:
        % vin: the input vector you want to map, range [min(vin),max(vin)]
        % rout: the range of the resulting vector
        % output:
        % vout: the resulting vector in range rout
        % usage:
        % >> v1 = linspace(-2,9,100);
        % >> v2 = linmap(v1,[-5,5]);
        %
        a = min(vin);
        b = max(vin);
        c = rout(1);
        d = rout(2);

        if (a == b) && (c == d)
            vout = vin;
        else
            vout = ((c+d) + (d-c)*((2*vin - (a+b))/(b-a)))/2;
        end

    end
end
% *** vgrid ***
%
% Val Schmidt
% Center for Coastal and Ocean Mapping
% University of New Hampshire
% Copyright 2010-2020, All rights reserved.

classdef vgrid

% TO DO:
%  Split properties into those the user can change and those he/she cannot.
%  Allow user to specify grid coordinates explicity.
properties

    cs        % cell size
    cinf      % cell influence
    type      % grid type

    xx        % x (Easting) coordinates of grid
    yy        % y (Northing) coordinates of grid
    nn        % Number of points contributing to grid cell
    zw        % Sum of product of gridded values and their weights for the cell
    ww        % Sum of weights for the cell.
    std       % Standard deviation of points in a cell (unweighted)
    
    Z         % Sequential estimator of depth for scalar (CUBE) or platlet methods.
    CZ        % Uncertainty or Covariance of Sequential Estimator.


end

methods

    function g = vgrid(cs,cinf,type,varargin)
        %% Function to create a new grid object. 
        %
        % INPUT:
        %   cs:    cell size (distance btwn nodes)
        %   cinf:  radial distance within which measurements contribute to
        %          a node.
        %   type:  grid type (see below)
        %   [OPTIONAL x,y and z must all be specified. w may or may not accompany them.]
        %   x:     x coordinates of input data.
        %   y:     y coordiantes of input data.
        %   z:     z coordinates of input data.
        %   w:     weights to be applied to each z value for weighted mean
        %          or parameter for other grid types i.e. Ngreaterthan w.
        %
        % Note:
        %   x,y and z must be the same size. w may be either unspecified or
        %   a scalar value (in either case all z values are equally
        %   weighted), or must be the same size as x, y, and z.
        %
        % Grid types:
        % mean:
        %   Average of the values. When w ~= 1, the mean is calculated by
        %   multipying each value in the cell by its weight divided by the sum of
        %   the weights in that cell.
        %
        % median:
        %   Calculates the median value for each grid cell.
        %
        % mode:
        %   Calculates the mode of the values for each grid cell.
        %
        % shoalest:
        %   Calculates the minimum value for each grid cell.
        %
        % deepest:
        %   Calculates the maximum value for each grid cell.
        %
        % stddev:
        %   Calculates the standard deviation of the values in each grid cell.
        %
        % stderr:
        %   Calculates the standard error of the values in each grid cell
        %   (stddev/N, where stddev is the standard deviation and N is the number
        %   of points falling in the cell)
        %
        % dwm:
        %   Calculates the distance weighted mean where each value in the cell is
        %   inversely weighted by the square if it's distance to the cell node.
        %
        % Nlowerthan:
        %   Calculates the number of points in the grid cell lower than some value,
        %   w.
        %
        % Ngreaterthan:
        %   Calculates the number of points greater than some value w.
        %
        % scalar_est:
        %   Scalar sequential estimator of depth a la CUBE. (Not yet
        %   implemented)
        %
        % platlet_est:
        %   Platlet sequential estimator of depth. (Not yet implemented)
        %
        % To Do:
        % - Rewrite mean function as a matrix operation to simplify the propagation
        % of uncertainty calcualtion. Actually this might be make more general such
        % that your pass a list of values, their uncertainty and weighting factors
        % and get back a mean and propagated uncertainty. This would allow
        % relatively simple incorporation of things like range weighting, footprint
        % weighting, gaussian weighting, etc.
        % - Add uncertainty to z input and propagate these through the
        % calculations.
        % - Add uncertainty to x and y inputs and propagate these through the
        % calculations (more difficult)
        % Rewrite a C mex function.
        %
        % Val Schmidt
        % CCOM/JHC
        % Copyright 2010-2020, All rights reserved.

        g.cs = cs;
        g.cinf = cinf;
        g.type = type;

        if ~isempty(varargin)
            g = g.add(varargin{:});
        end

    end

    g = add(g,x,y,z,w);       % Add data to the grid
    z = zz(g);                % Calculate the z values for the grid.
    g = rotate(g,r,p,h);      % Rotate a grid.


end


end
function g = add(g,x,y,z,w,varargin)
%% An incremental gridding function
% Arguments:
% x:   x-coordinates
% y:   y-coordiantes
% z:   z-scalar values to grid
% w:   w-weight applied to each point (size of x or 1 for no weighting)
%      When 'type' = Nlowerthan or Ngreaterthan, w is the threshold value
%      When 'type' = distance weighted mean, distance = R^w
% cs:  grid cell size
% cinf: cell influence 
% type: type of grid (see below)
% g:   Optional previously defined grid structure
%
% Output:
% g.xx: vector of grid cell x coordinates.
% g.yy: vector of grid cell y coordiantes.
% g.zz: 2D matrix of grided values times their weights.
% g.nn: 2D matrix containing the number of points in each grid cell.
% g.ww: sum of weights of items in the grid cell
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
% Center for Coastal and Ocean Mapping
% University of New Hampshire
% Copyright 2010-2020, All rights reserved.

debug = 0;
newgrid = 0;

x = x(:);
y = y(:);
z = z(:);
w = w(:);

% Weight cannot be zero. 
if ~isscalar(w)
    if sum(w == 0)
        warning('Found zero weights. Weights cannot be 0. Setting to 1e-20')
        w(w==0) = 1e-20;
    end
end

% Format for new grid
if isempty(g.zw)
    
    newgrid = 1;
   
    
    % Get coordinates of grid cells
    g.xx = min(x(:)):g.cs:(max(x(:))+g.cs);
    g.yy = min(y(:)):g.cs:(max(y(:))+g.cs);
    
    % Ask if the grid is REALLY huge.
    if ~(gridsizesanitycheck(g.xx) | gridsizesanitycheck(g.yy))
        return
    end
    
    

else
   
    sprintf('NOTE: The algorithm for the standard deviation of gridding iterations\n');
    sprintf('approximates the standard deviation by the mean of the square of the\n');
    sprintf('standard deviation of the last iteration and that of the\n');
    sprintf('current iteration.\n');
    
    
    % Extend the current grid.
    minx = min(x(:));
    miny = min(y(:));
    maxx = max(x(:));
    maxy = max(y(:));
    
    if minx < g.xx(1)
        dx = minx:g.cs:(g.xx(1)-g.cs);
        if ~gridsizesanitycheck(g.xx)
            return
        end
        g.xx = [dx g.xx];
        g.zw = [nan(length(g.yy),length(dx)) g.zw];
        g.nn = [nan(length(g.yy),length(dx)) g.nn];
        g.ww = [nan(length(g.yy),length(dx)) g.ww];
        g.std = [nan(length(g.yy),length(dx)) g.std];
        if strcmp(g.type,'platelest_est')
            
            % A tricky way to create a cell array length(g.yy) x
            % length(g.xx) filled with 3x1 vectors of nans. 
            tmp = mat2cell(nan(length(g.yy)*3,length(dx)),...
                ones(length(g.yy),1)*3,ones(length(dx),1));            
            g.Z = {tmp g.Z};
            % Same to init 3x3 nan cell array.
            tmp = mat2cell(nan(length(g.yy)*3,length(dx)*3),...
                ones(length(g.yy),1)*3,ones(length(dx),1)*3);

            g.ZC = {tmp g.ZC};
            
        end

    end

    if maxx > g.xx(end)
        dx = (g.xx(end)+g.cs):g.cs:maxx;
        if ~gridsizesanitycheck(g.xx)
            return
        end
        g.xx = [g.xx dx];
        g.zw = [g.zw nan(length(g.yy),length(dx))];
        g.nn = [g.nn nan(length(g.yy),length(dx))];
        g.ww = [g.ww nan(length(g.yy),length(dx))];
        g.std = [g.std nan(length(g.yy),length(dx))];
        
        if strcmp(g.type,'platelest_est')
            
            % A tricky way to create a cell array length(g.yy) x
            % length(g.xx) filled with 3x1 vectors of nans. 
            tmp = mat2cell(nan(length(g.yy)*3,length(dx)),...
                ones(length(g.yy),1)*3,ones(length(dx),1));            
            g.Z = {g.Z tmp};
            % Same to init 3x3 nan cell array.
            tmp = mat2cell(nan(length(g.yy)*3,length(dx)*3),...
                ones(length(g.yy),1)*3,ones(length(dx),1)*3);

            g.ZC = {g.ZC tmp};
            
        end
    end
    if miny < g.yy(1)
        dy = miny:g.cs:(g.yy(1)-g.cs);
        if ~gridsizesanitycheck(g.yy)
            return
        end
        g.yy = [dy g.yy];
        g.zw = [nan(length(dy), length(g.xx)); g.zw];
        g.nn = [nan(length(dy), length(g.xx)); g.nn];
        g.ww = [nan(length(dy), length(g.xx)); g.ww];
        g.std = [nan(length(dy), length(g.xx)); g.std];
        
        if strcmp(g.type,'platelest_est')
            
            % A tricky way to create a cell array length(g.yy) x
            % length(g.xx) filled with 3x1 vectors of nans. 
            tmp = mat2cell(nan(length(dy)*3,length(g.xx)),...
                ones(length(dy),1)*3,ones(length(g.xx),1));            
            g.Z = {tmp; g.Z};
            % Same to init 3x3 nan cell array.
            tmp = mat2cell(nan(length(dy)*3,length(g.xx)*3),...
                ones(length(dy),1)*3,ones(length(g.xx),1)*3);

            g.ZC = {tmp; g.ZC};
            
        end
    end
    if maxy > g.yy(end)
        dy = (g.yy(end)+g.cs):g.cs:maxy;
        if ~gridsizesanitycheck(g.yy)
            return
        end
        g.yy = [g.yy dy];
        g.zw = [g.zw; nan(length(dy), length(g.xx))];
        g.nn = [g.nn; nan(length(dy), length(g.xx))];
        g.ww = [g.ww; nan(length(dy), length(g.xx))];
        g.std = [g.std; nan(length(dy), length(g.xx))];
        
        if strcmp(g.type,'platelest_est')
            
            % A tricky way to create a cell array length(g.yy) x
            % length(g.xx) filled with 3x1 vectors of nans. 
            tmp = mat2cell(nan(length(dy)*3,length(g.xx)),...
                ones(length(dy),1)*3,ones(length(g.xx),1));            
            g.Z = {g.Z; tmp};
            % Same to init 3x3 nan cell array.
            tmp = mat2cell(nan(length(dy)*3,length(g.xx)*3),...
                ones(length(dy),1)*3,ones(length(g.xx),1)*3);

            g.ZC = {g.ZC; tmp};
            
        end
    end
    
    
end
    
% Not using this at the moment. Could return XX and YY instead of xx and
% yy, but it takes more memory to do so. 
%[XX YY] = meshgrid(xx,yy);

% Get boundaries of the grid cells
% xxlb = g.xx-g.cs/2;
% xxup = g.xx+g.cs/2;
% yylb = g.yy-g.cs/2;
% yyup = g.yy+g.cs/2;

xxlb = g.xx-g.cinf;
xxup = g.xx+g.cinf;
yylb = g.yy-g.cinf;
yyup = g.yy+g.cinf;

% Initialize grid matrix
grows = length(g.yy);
gcols = length(g.xx);

%zz = cl_memmapfile('Value',0,'Size',[grows gcols]);
%nn = cl_memmapfile('Value',0,'Size',[grows gcols]);
doindices = 0;

if newgrid
    g.zw = nan(grows,gcols);
    g.nn = nan(grows,gcols);
    g.ww = nan(grows,gcols);
    g.std = nan(grows,gcols);
    
    % Initialize the values for the platelet estimation. 
    if strcmp(g.type,'platlet_est')
        g.Z = mat2cell(nan(3*grows,gcols),ones(grows,1)*3,ones(gols,1));
        g.CZ = mat2cell(nan(3*grows,3*gcols),ones(grows,1)*3,ones(gols,1)*3);
    end

    % This is for experimental returning of indices in which individual
    % data points contriubted. 
%     if strcmp(g.type,'indices')
%         g.type = 'mean'
%         nn = nan(size(z));
%         doindices = 1;
%     else
%         nn = nan(grows,gcols);
%     end
    
end


if debug
    [xxx yyy] = meshgrid(g.xx,g.yy);
    plot(xxx(:),yyy(:),'.k')
    hold on
    plot(x,y,'.');
    hold off
end

cinf2 = g.cinf.^2;
for idx = 1:grows
    
    % We need to search through all the data efficiently to determine
    % indices for points that will contribute to a grid node. Those that
    % contribute are ones that fall within the "cell influence" (cinf).
    % Thse are the ones that meet the criteria: 
    %
    % sqrt( (x-xo).^2 + (y-yo).^2 ) < cinf
    %
    % Reorganizing....
    %
    % (x-xo)^2 + (y-yo)^2 < cinf^2
    %
    % This will never be true when either term of the lhs is >= cinf^2. So
    % we reduce the search by doing these piece-meal. 
    
    
    % Here we find the y data values within cinf of the grid node.
    %yidx = find ( (y-yy(idx)).^2 < cinf2);
    % Same find statement reworked to support cl_memmapfile.
    ddy = (y-g.yy(idx)).^2;
    %yidx = find ( y < g.yy(idx) + g.cinf & y > g.yy(idx) - g.cinf);
    yidx = find( ddy < cinf2);
    
    if isempty(yidx)
        continue
    end
    
    for jdx = 1:gcols
        % Then go through each cell of that row, and look for x-values that
        % also are in the cell.
        
        %xidx = find( x(yidx) > xxlb(jdx) & x(yidx) < xxup(jdx) );
        %xidx = find( (x(yidx) - xx(jdx)).^2 < cinf2);
        % Same find statement reworked to support cl_memmapfile.
        %xidx = find( x(yidx) < g.xx(jdx) + g.cinf & x(yidx) > g.xx(jdx) - g.cinf);
        xidx = find( (x(yidx) - g.xx(jdx)).^2 < cinf2 - ddy(yidx));
        
        if isempty(xidx)
            continue
        end
        if debug
            hold on
            plot(xxx(idx,jdx),yyy(idx,jdx),'om','linewidth',2)
            plot(x(yidx(xidx)),y(yidx(xidx)),'.g');
            hold off
            keyboard
        end
        
        if strcmp(g.type,'dwm')
            % Calculate distance between points for distance-weighted-mean
            % algorithm. In this case W is the exponent
            R = ((xx(jdx)-x(yidx(xidx))).^2 + ( yy(jdx)-y(yidx(xidx))).^2).^(w/2);
        end
        
        %drawnow;
%        g.nn(idx,jdx) = nansum([g.nn(idx,jdx) length(yidx(xidx))]);
        
        % If doindices is set, nn returns the index of the grid node to
        % which the point was last attributed. Otherwise it returns the
        % number of data points attributed to the grid node.
        if ~doindices
            g.nn(idx,jdx) = nansum([g.nn(idx,jdx) length(xidx)]);
        else
            g.nn(yidx(xidx)) = idx*(colz-1)+jdx;          
        end
        
        switch g.type
            case 'mean'
                if isscalar(w) % no weighting.
                    g.zw(idx,jdx) = nansum([g.zw(idx,jdx); z(yidx(xidx))]);
                    g.ww(idx,jdx) = g.nn(idx,jdx);  % Shouldn't this be multiplied by w???
                    g.std(idx,jdx) = sqrt(nanmean([g.std(idx,jdx); nanstd(z(yidx(xidx)))].^2));
                else
                    % weighted gridding = value times the weight divided by the sum of the weights
                    g.zw(idx,jdx) = nansum([g.zw(idx,jdx); z(yidx(xidx)).*w(yidx(xidx))]);
                    g.ww(idx,jdx) = nansum([g.ww(idx,jdx); w(yidx(xidx))]);
                    g.std(idx,jdx) = sqrt(nanmean([g.std(idx,jdx); nanstd(z(yidx(xidx)))].^2));
                end
            case 'meanwithoutlierrejection'
                if isscalar(w) % no weighting.
                   error('w must be a vector of uncertainty values for these soundings')
                else
                    % This is going to be odd and may produce odd results.
                    % We are going to assume that, as points are added,
                    % there are enough within each batch to provide some
                    % statisical power. So if the cell is empty, we'll
                    % calculate the mean and stddev, omit points greater
                    % than x-sigma, then calculate the mean again. If the
                    % cell is not empty, we're going to assume the compare
                    % the points to the current depth and uncertainty value
                    % and omit the ones greater than x-sigma. 
                    if isnan(z(idx,jdx))
                        z(idx,jdx) = z(yidx(xidx(1)));
                    end
                    
                    for ii=1:length(xidx)
                       
                        
                        
                        
                    end
                        
                        
                    
                    
                    tmpzw = nansum([g.zw(idx,jdx); z(yidx(xidx)).*w(yidx(xidx))]);
                    tmpww = nansum([g.ww(idx,jdx); w(yidx(xidx))]);
                    tmpstd = sqrt(nanmean([g.std(idx,jdx); nanstd(z(yidx(xidx)))].^2));
                    
                    % Loop for outliers.
                    M = abs(tmpzw./tmpww - z(yidx(xidx))) < tmpstd*3;
                    % If we find them, adjust our calculation omittin them.
                    if any(~M)
                        itmp = yidx(xidx);
                        g.zw(idx,jdx) = nansum([g.zw(idx,jdx); z(itmp(M)).*w(itmp(M))]);
                        g.ww(idx,jdx) = nansum([g.ww(idx,jdx); w(itmp(M))]);
                        g.std(idx,jdx) = sqrt(nanmean([g.std(idx,jdx); nanstd(z(itmp(M)))].^2));
                    else
                        % otherwise accepth it.
                        g.zw(idx,jdx) = tmpzw;
                        g.ww(idx,jdx) = tmpww;
                        g.std(idx,jdx) = tmpstd;
                    end
                        
                end
                
            case 'median'
                % Not sure this makes sense, but here I'm calculating the
                % current median (or taking what was calculated last) and
                % the number of points between the max and the median. Then
                % I 
%                 g.max(idx,jdx) = nanmax(z(yidx(xidx)));
%                 g.med(idx,jdx) = nanmedian(z(yidx(xidx)));
%                 g.m_idx(idx,jdx) = find(g.med(idx,jdx)==z(yidx(xidx)),1,'first');
%                 
                g.zw(idx,jdx) = nanmedian([g.zw(idx,jdx) nanmedian(z(yidx(xidx)))]);  
                % g.zz(idx,jdx) = nanmedian(z(yidx(xidx)));
                g.ww(idx,jdx) = 1;
                % This doesn't make sense for incremental uncertainty, but
                % it'll be a reasonable first go.
                g.std(idx,jdx) = sqrt(nanmean([g.std(idx,jdx); nanstd(z(yidx(xidx)))].^2));

            case 'mode'
                warning('Incremental mode not implmented')
                zz(idx,jdx) = mode(z(yidx(xidx)));
                g.ww(idx,jdx) = 1;
            case 'shoalest'
                g.zw(idx,jdx) = nanmin([g.zw(idx,jdx); z(yidx(xidx))]);
                g.ww(idx,jdx) = 1;
            case 'deepest'
                g.zw(idx,jdx) = nanmax([g.zw(idx,jdx); z(yidx(xidx))]);
                g.ww(idx,jdx) = 1;
            case 'stddev'
                % If I do this with weights I'll have implemented CUBE!
                 g.std(idx,jdx) = sqrt(nanmean([g.zw(idx,jdx).^2; nanstd(z(yidx(xidx)).^2)]));
                %g.std(idx,jdx) = nanstd(z(yidx(xidx)));
                %g.ww(idx,jdx) = 1;
            case 'stderr'
                g.zw(idx,jdx) = nanstd(z(yidx(xidx)))./sqrt(g.nn(idx,jdx));
                g.ww(idx,jdx) = 1;
            case 'dwm'
                % experimental
                g.zw(idx,jdx) = g.zw(idx,jdx) + nansum( z(yidx(xidx))).*(1./R);
                g.ww(idx,jdx) = g.ww(idx,jdx) + sum(1./R);
            case 'Nlowerthan'
                g.zw(idx,jdx) = sum( z(yidx(xidx)) < g.ww ) + g.zw(idx,jdx);
            case 'Ngreaterthan'
                g.zw(idx,jdx) = sum( z(yidx(xidx)) > g.ww ) + g.zw(idx,jdx);
                g.ww(idx,jdx) = 1;
                
            case 'scalar_est'
                error('Grid type not yet supported.')
            case 'platlet_est'
                error('Grid type not yet supported.')
                
                % If no value is yet assigned and N <3, average the two
                % values.
                if g.Z(idx,jdx) == nan(3,1) && numel(z(yidx(xidx))) < 3
                   % FIX: Should this be the weighted sum instead.
                   g.Z(idx,jdx) = [0 0 nansum(z(yidx(xidx)).*w(yidx(xidx)))./sum(w(yidx(xidx)))];
                   % FIX: This is not the correct way to calculate the
                   % uncertainty of a weighted sum of values. 
                   g.CZ{idx,jdx} = eye(3).*std(z(yidx(xidx)));
                end
                
                % If no value is yet assigned and N >= 3, calculate a least
                % squares solution.
                
                % Create a measurement design matrix, columes of x's, y's
                % and 1's. Note we subtract off the location of the grid
                % node center.
                h = [x(yidx(xidx)) y(yidx(xidx)) ones(numel(yidx(xidx)),1)] ...
                    - [g.xx(idx,jdx) g.yy(idx,jdx) 0]
                C = diag(w(yidx(xidx)));
                
                Z{idx,jdx} = inv(h'*C*h)*h*C*z(yidx(xidx));
                C{idx,jdx} = inv(h'*C*h);  % FIX: Check this!
                
                % If a value is assigned already, step through each new
                % measurement and update the value with a sequential linear
                % estimator. 
                
                
                
            otherwise
                error('Grid type not supported.')
                return;
        end
    end
    if debug
        pause(1)
    end
end

function Q = gridsizesanitycheck(mm)
%% Check to see if the grid size is going to be REALLY large and warn the user.
Q = 1;
if length(mm) > 1e4
    disp(['Grid bounds are: ' num2str(min(x)) ':' num2str(max(x)) ',' ...
        num2str(min(y)) ':' num2str(max(y)) ' Grid size is ' num2str(length(g.xx)) ...
        ':' num2str(length(g.yy))])
    R = input('Do you want to continue? (y)/n  >','s');
    if R == 'n'
        Q = 0;   
    end
end

end


end






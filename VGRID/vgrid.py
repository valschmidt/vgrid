#!/usr/bin/env python
#
# ---- VGRID ----
# A Class for Incremental Gridding
#
# Val Schmidt
# Center for Coastal and Ocean Mapping
# Univeristy of New Hampshire
# Copyright, 2018-2019
# All Rights Reserved


import numpy as np
from line_profiler import LineProfiler
import sys

class vgrid():
    ''' A class for gridding of x,y,z data.
    
    See vgrid.add() for details on usage.
    '''

    def __init__(self, cs=1.0, cinf=1.0, type='mean'):

        self.cs = cs        # cell size
        self.cinf = cinf    # cell influence
        self.type = type    # grid type

        self.xx = None           # x (Easting) coordinates of grid
        self.yy = None           # y (Nothing) coordinates of grid
        self.zw = None           # Sum of the product of the gridded values and their weights for the grid cell.
        self.ww = None           # Sum of weights
        self.nn = None           # Number of points contributing to grid cell

        self.varw = None          # Sum of the square of the difference of the gridded values and the estimated mean, times their weights.
        self.Z = None            # Sequential estimator of depth for scalar (CUBE) or platlet methods.
        self.CZ = None
        
        ### Utility variables used internally. ###
        # New values to incorporate into grid.
        self._x = None
        self._y = None
        self._z = None
        self._w = None
        self._II = None      # Indices of values to add for node under consideration.
        # (These are not used as class variables currently.)
        self._idx = None     # row grid node indiex for node under consideration.
        self._jdx = None     # column grid node indiex for node under consideation.

    def zz(self):
        ''' Calculate the z values for the grid.'''
        return self.zw / self.ww

    def mean(self,idx,jdx):
        '''Mean gridding algorithm.

	    vgrid implemnets incremental gridding where possible. 
	    To do this, the sum of the product of the weights and z values are retained
        in addition to the sum of the weights. Then method zz() calculates the 
        quotient of the two to obtain the actual weighted mean z values. Note that
        when all weights are one, (or if w is set to 1 for shorthand), a standard
        mean is calculated.

        Variance is calcualted in a similar way. In this case the sum of w*(z_i - mu)^2
        is calculated and stored for each grid node, where z_i is the value to be gridded
        and mu is the mean of the grid node calculated thus far.  Then this sum 
        is divided by the sum of the weights to get the final estimated variance. As the
        mean of the grid node approaches the true mean, this value should approach the 
        true variance. 
        '''

        # Non-weighted gridding.
        if self._w.size == 1:
            self.zw[idx, jdx] = np.nansum(np.concatenate((self._z[self._II], [self.zw[idx, jdx]])))
            self.ww[idx, jdx] = self.nn[idx, jdx]
            self.varw[idx, jdx] = np.nansum(np.concatenate((np.power( (self._z[self._II] - self.zw[idx,jdx]/self.nn[idx,jdx]), 2),[self.varw[idx, jdx]])))
        else:
            # Weighted gridding. Sum of value times the weight divided by the 
            # sum of the weights.
            # The strategy taken here is to retain the sum of the values times the weights, and also
            # the sum of the weights. Then when the weighted mean is requested the calling function 
            # divides the these two values. This strategy allows incremental addition of data to the grid.
            #
            # The coding strategy below is to append the new points to the existing point in a list
            # and then call nansum to add them up. 
            #
            # Q: Note: A dot-product might be quicker, but there is no dot-product that will produce a
            # non-nan result if one of the values is nan, which is desired here.
            self.zw[idx, jdx] = np.nansum(np.append(self.zw[idx, jdx], self._z[self._II] * self._w[self._II]))
            self.ww[idx, jdx] = np.nansum(np.append(self.ww[idx, jdx], self._w[self._II]))
            self.varw[idx, jdx] = np.nansum(np.append(np.power( (self._z[self._II] - self.zw[idx,jdx]/self.ww[idx,jdx]),2)
                                                      , self.varw[idx, jdx] ))

    def var(self):
        ''' Calculate the variance'''
        return self.varw/self.ww

    def std(self):
        '''Calculate the standard deviation'''
        return np.sqrt(self.var())

    def meanwithoutlierrejection(self):
        ''' TO DO: Calculate the mean, rejecting values that exceed 3-sigma
        from existing estimate.'''
        pass

    def median(self, idx, jdx):
        ''' Calculate the median value in each grid cell.
        
        The method used here to provide a "running median" is for each add(),
        calculate the average of the existing value with the median of the
        new points. This method works reasonably well, but can produce
        inferior results if a single add() contains only outliers and their
        are insufficient additional adds to constrain it.'''
        self.zw[idx, jdx] = np.nanmean(
            np.append(self.zw[idx, jdx], np.nanmedian(self._z[self._II])))
        self.ww[idx, jdx] = 1
        self.varw[idx, jdx] = np.nansum(np.append(
            np.power((self._z[self._II] - self.zw[idx, jdx]/self.ww[idx, jdx]),
                     2),
            self.varw[idx, jdx]))
        pass

    def gridsizesanitycheck(self, M):
        '''Check to see if the grid size is going to be REALLY large. '''
        if M.__len__() > 1e4:
            return False
        else:
            return True
        
    def create_new_grid(self):
        ''' Create a new empty grid.'''
        
        self.xx = np.arange(min(self._x), max(self._x)+self.cs, self.cs)
        self.yy = np.arange(min(self._y), max(self._y)+self.cs, self.cs)

        if not (self.gridsizesanitycheck(self.xx) and
                self.gridsizesanitycheck(self.yy)):
            print('Grid size is too large.')
            return
        
        # Initialize grid.
        self.zw = np.empty((self.yy.size, self.xx.size))
        self.zw.fill(np.nan)
        self.nn = np.copy(self.zw)
        self.ww = np.copy(self.zw)
        self.varw = np.copy(self.zw)
        
    def expand_grid(self):
        minx = min(self._x)
        miny = min(self._y)
        maxx = max(self._x)
        maxy = max(self._y)

        if minx < self.xx[0]:
            dx = np.arange(minx, self.xx[0] - self.cs, self.cs)
            self.xx = np.concatenate((dx,self.xx))
            # Create new space
            tmp = np.empty((self.yy.size,dx.size))
            tmp.fill(np.nan)
            # Tack it on.
            self.zw = np.concatenate((np.copy(tmp), self.zw), axis=1)
            self.nn = np.concatenate((np.copy(tmp), self.nn), axis=1)
            self.ww = np.concatenate((np.copy(tmp), self.ww), axis=1)
            self.varw = np.concatenate((np.copy(tmp), self.varw), axis=1)

        # FIX: Support depth/platelet estimates here, tbd

        if maxx > self.xx[-1]:
            dx = np.arange(self.xx[-1]+self.cs,maxx,self.cs)
            self.xx = np.concatenate((self.xx,dx))
            # Create new space
            tmp = np.empty((self.yy.size,dx.size))
            tmp.fill(np.nan)
            # Tack it on.
            self.zw = np.concatenate((self.zw,np.copy(tmp)),axis=1)
            self.nn = np.concatenate((self.nn,np.copy(tmp)),axis=1)
            self.ww = np.concatenate((self.ww,np.copy(tmp)),axis=1)
            self.varw = np.concatenate((self.varw,np.copy(tmp)),axis=1)

        if miny < self.yy[0]:
            dy = np.arange(miny,self.yy[0]-self.cs,self.cs)
            self.yy = np.concatenate((dy,self.yy))
            tmp = np.empty((dy.size,self.xx.size))
            tmp.fill(np.nan)
            self.zw = np.concatenate((np.copy(tmp), self.zw),axis=0)
            self.nn = np.concatenate((np.copy(tmp), self.nn),axis=0)
            self.ww = np.concatenate((np.copy(tmp), self.ww),axis=0)
            self.varw = np.concatenate((np.copy(tmp), self.varw),axis=0)

        if maxy > self.yy[-1]:
            dy = np.arange(self.yy[-1]+self.cs,maxy,self.cs)
            self.yy = np.concatenate((self.yy,dy))
            tmp = np.empty((dy.size,self.xx.size))
            tmp.fill(np.nan)
            self.zw = np.concatenate((self.zw,np.copy(tmp)), axis=0)
            self.nn = np.concatenate((self.nn,np.copy(tmp)), axis=0)
            self.ww = np.concatenate((self.ww,np.copy(tmp)), axis=0)
            self.varw = np.concatenate((self.varw,np.copy(tmp)), axis=0)


    def add(self, x, y, z, w):
        ''' An incremental gridding function

        Arguments:
        x:   x-coordinates
        y:   y-coordiantes
        z:   z-scalar values to grid
        w:   w-weight applied to each point (size of x or 1 for no weighting)
             When 'type' = Nlowerthan or Ngreaterthan, w is the threshold value
             When 'type' = distance weighted mean, distance = R^w
        cs:  grid cell size
        cinf: cell influence
        type: type of grid (see below)

        Output:
        g.xx: vector of grid cell x coordinates.
        g.yy: vector of grid cell y coordiantes.
        g.zz: 2D matrix of grided values times their weights.
        g.nn: 2D matrix containing the number of points in each grid cell.
        g.ww: sum of weights of items in the grid cell

        %
        % Grid types:
        % mean:
        %   Average of the values. When w != 1, the mean is calculated by
        %   multipying each value in the cell by its weight divided by the sum
        %   of the weights in that cell.
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
        % CCOM/JHC
        % 2018, 2019
        '''

        # Force everything to match.
        if np.isscalar(x) or np.isscalar(y) or np.isscalar(z):
            print('X, Y, or Z is scalar - must be numpy array.')
            sys.exit()

        self._x = x.ravel()
        self._y = y.ravel()
        self._z = z.ravel()
        if not np.isscalar(w):
            self._w = w.ravel()
        else:
            self._w = np.array(w)

        # Weight cannot be zero.
        if self._w.size != 1:
            if sum(self._w == 0):
                print("Found zero weights. Weights cannot be zero.")
                print("Setting to 1e-20.")
                self._w[self._w == 0] = 1e-20

        # Set up new grid, or extend the existing grid if necessary.
        if self.zw is None:
            self.create_new_grid()

        else:
            self.expand_grid()

        grows = self.yy.size
        gcols = self.xx.size

        doindices = 0

        cinf2 = self.cinf**2

        # Go through the rows of the grid..
        for idx in range(grows):
            '''
            We need to search through all the data efficiently to determine
            indices for points that will contribute to a grid node. Those that
            contribute are ones that fall within the "cell influence" (cinf).
            Thse are the ones that meet the criteria:

            sqrt( (x-xo).^2 + (y-yo).^2 ) < cinf

            Squaring both sides....

            (x-xo)^2 + (y-yo)^2 < cinf^2

            This will never be true when either term of the lhs is >= cinf^2.
            So we reduce the search by doing these piece-meal. '''



            # Here we find the y data values within cinf of the grid node
            ddy = (self._y - self.yy[idx])**2
            #yidx = np.flatnonzero(ddy < cinf2)
            yidx = np.flatnonzero(ddy < cinf2)
            #yidx = ddy < cinf2

            # If there are none, then don't bother with further calculations.
            if yidx.size == 0:
                continue

            # Then go through each cell of that row, and look for x - values that also are in the cell.
            # But first pre-calculate a vector of terms that will be needed for every evaluation.
            xtest = cinf2 - ddy[yidx]
            for jdx in range(gcols):
                xidx = np.flatnonzero( (self._x[yidx] - self.xx[jdx])**2 < xtest )

                # If there are none of these then there is nothing to do to the grid node.
                if xidx.size == 0:
                    continue
                
                # Set the indices of the values to be add to this grid node.
                self._II = yidx[xidx]

                if self.type == 'dwm':
                    # Calculate distance between points and grid node for distance-weighted mean.
                    # In the case, w is the exponent.
                    R = ((self.xx[jdx] - self._x(self.II))**2 +
                         (self.yy[jdx]-self._y[self.II])**2)**(self._w/2.0)

                if not doindices:
                    self.nn[idx,jdx] = np.nansum(np.append(self.nn[idx,jdx], xidx.size))
                else:
                    self.nn[idx,jdx] = idx*(gcols-1) + jdx

                #print('INDEXES: %s' % ','.join(map(str,yidx[xidx])))
                #print('VALUES: %s' % ','.join(map(str,z[yidx[xidx]])))

                if self.type == 'mean':
                    self.mean(idx,jdx)
                    
                if self.type == "median":
                    self.median(idx,jdx)



        def rotate(self):
            pass

    def pcolor(self,*kwargs):

        from matplotlib import pyplot as plt

        plt.pcolor(self.xx,self.yy,self.zz(),*kwargs)
        #plt.colorbar()
        plt.ion()
        plt.show()
        #plt.draw()
        plt.pause(0.001)



if __name__=='__main__':

    profileON = True

    def gridTest(N = 2, ProfileON = False):
        ''' Method to test gridding.'''

        print("N=%d" % N)
        # Generate data. 
        x = np.random.random((N,1))*100
        y = np.random.random((N,1))*100
        z = np.exp( np.sqrt((x-50.)**2 + (y-50.)**2)/50)

        # Generate grid.
        G = vgrid(1,1,'mean')
        if profileON:
            print("Profiling on.")
            lp = LineProfiler()
            GAddProfiled = lp(G.add)
            lp.add_function(G.mean)
            GAddProfiled(x,y,z,1)
            return (G,lp)
        else:
            G.add(x,y,z,1)
            return G
  
    ##  Gridding test script:

    for N in [1000, 5000, 10000, 20000, 50000, 100000, 1000000]:
        if profileON:
            GG,LP = gridTest(N,True)
            LP.print_stats()
        else:
            GG = gridTest(N,False) 

        # Plot test.
        GG.pcolor()




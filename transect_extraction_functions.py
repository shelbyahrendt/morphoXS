from shapely.geometry import LineString, Point
import numpy as np

def addpts(linestring, dl):
    '''
    This function generate points at fixed interval dl in line units along a shapely Linestring

    Parameters
    ------------
    linestring: a Shapely LineString geometry
    dl: distance interval to generate points

    Output
    --------
    l: list of distance along line
    xpts: list of x coordinates
    ypts: list of y coordinates

    Note: this function was inspired by David Shean's pygeotools function 'line2pts' https://github.com/dshean/pygeotools/tree/master
    '''

    #define lists to hold output 'x,y' points and distance 'l' along line
    #populate with first x,y verticies and '0' for first distance
    l = [0]
    xpts = [linestring.coords[0][0]]
    ypts = [linestring.coords[0][1]]
    #remainder variable to carry little distances between linestring points
    rem_l = 0

    #this iterates through existing verticies in line
    for ((x1,y1), (x2,y2)) in zip(linestring.coords[:-1], linestring.coords[1:]):
        #get distance between existing verticies
        dL = ((x2-x1)**2 + (y2-y1)**2)**(1/2)
        #get how many points 'dl' apart we can put here
        n = int((rem_l + dL) / dl)
        
        #loop through points and add them if n>0
        for pt in range(n):
            if pt == 0:
                #apply remainder dx and dy bits first
                rem_dx = (x2-x1) * (dl-rem_l) / dL
                rem_dy = (y2-y1) * (dl-rem_l) / dL
                xpts += [x1 + rem_dx]
                ypts += [y1 + rem_dy]
                l += [l[-1] + dl] #insert an 'l' value since we added a point
            else:
                #convert dl to dx for segment
                dx = (x2-x1) * dl/dL
                dy = (y2-y1) * dl/dL
                xpts.append(xpts[-1] + dx)
                ypts.append(ypts[-1] + dy)
                l += [l[-1] + dl] #insert an 'l' value since we added a point
            
        # update remainder between existing LineString verticies
        rem_l += dL - n*dl
    return(l, xpts, ypts)

def make_perp_transects(xpts, ypts, L, DEM_res):
    '''
    This function creates perpindicular transects of length L
    at certain locations along a channel centerline defined by xpts y pts
    (Pair with addpts function to get equi-spaced points along centerline)
    
    Parameters
    ------------
    xpts: list of points defining center of transect x locations (compatible with 'addpts' output)
    ypts: list of points defining center of transect y locations (compatible with 'addpts' output)
    L: length of transect
    DEM_res: resolution to extract topography at each transect(set to DEM res for finest interpolation)

    Output
    --------
    transects: a dictionary defining the points along each transect of form ([x], [y], [d]) where
    x is a list of x locations along each transect in native coordinate system units
    y is a list of y locations along each transect in native coordinate system units
    d is a list of the distance along each transect
    This dictionary is keyed as transects[i] where 'i' is the index of a transect starting at '1'
    
    '''
    
    
    # transects = [] #list of tuples containing x,y arrays of points within each transect
    transects = {} #dictionary to hold transects
    n_vector = [] #normal vector defining the direction of each transect (important for unidirectional river flow)
    
    for i in range(1, len(xpts)-1): #n transects is two less than n points since using midpoint
            
        #calculate slope of transect for surrounding points
        S = (ypts[i+1] - ypts[i-1]) / (xpts[i+1] - xpts[i-1])
        St = -1/S
    
        #create distances to shift x and y pts away from line vertex
        x_shift = L / (2 * np.sqrt(1 + (St)**2)) 
        y_shift = St * x_shift
        
        #create endpoints of transect by moving points away from second point on centerline vector
        x1, y1 = (xpts[i] + (x_shift)), (ypts[i] + (y_shift)) 
        x2, y2 = (xpts[i] - (x_shift)), (ypts[i] - (y_shift))
        
        
        #create shapely line from this (uses shapely LineString and Point objects)
        # if statement keeps 0 point of each transect on left bank
        if ypts[i+1] > ypts[i-1]: #if river is flowing in positive y direction
            line = LineString([Point(x2, y2), Point(x1, y1)])
        elif ypts[i+1] < ypts[i-1]: #if river is flowing in negative y direction
            line = LineString([Point(x1, y1), Point(x2, y2)])
            
        #put intermediate points on line at resolution of DEM using line2pts function
        trans_l, trans_x, trans_y = addpts(line, DEM_res)
        # transects.append((trans_x, trans_y, trans_l))
        transects[i-1] = {'x': trans_x,
                         'y': trans_y,
                         'd': trans_l}
    
    return transects



def calc_curvature(xpts, ypts):
    '''
    This function computes curvature along a centerline given a list of x & y points
    
    Parameters
    ------------
    xpts: list of points defining center of transect x locations (compatible with 'addpts' output)
    ypts: list of points defining center of transect y locations (compatible with 'addpts' output)

    Output
    --------
    kappa: a list of curvature evaulated at x,y points in [1/L] units where 'L' is length in x,y units
    rc: radius of curvature in [L] units where 'L' is length in x,y units
    
    '''
    # first derivatives of x and y
    dx = np.gradient(xpts)
    dy = np.gradient(ypts)
    #second derivateives of x and y
    ddx = np.gradient(dx)
    ddy = np.gradient(dy)
    #curvature
    kappa = (dx*ddy - dy*ddx) / ((dx**2 + dy**2)**1.5)
    rc = 1/kappa

    return kappa, rc 
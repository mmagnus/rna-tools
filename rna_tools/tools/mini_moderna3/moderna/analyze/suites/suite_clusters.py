#!/usr/bin/python
#
# suite_clusterss.py
#
# Parses the set of torsion angles provided by Jane Richardson
#

__author__ = "Kristian Rother, Raphael Bauer"
__credits__ = ["Raphael Bauer","Markus Weber","Marcin Domagalski","Magdalena Musielak", "Janusz Bujnicki", "Marie Curie"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"


from math import pi,cos

EMPTY_ANGLES = [-1.0]*7

POWER = 3.0

NORMAL_WEIGHTS = [28.0, 60.0, 55.0, 50.0, 70.0, 35.0, 28.0]
SATELLITE_WEIGHTS = [28.0, 50.0, 50.0, 45.0, 60.0, 35.0, 28.0]

# adjusted radii for satellite-dominant comparisons.
# satnam:((1 s2 s3 4 s5 6 7),(1 d2 d3 4 d5 6 7 ,domnam))
DOM_SAT_EMPIRICAL = {
   "1m":((0.0, 0.0, 0.0,0.0,32.0,0.0,0.0),(0.0, 0.0, 0.0,0.0,64.0,0.0,0.0 ,"1a")),
   "1L":((0.0,18.0, 0.0,0.0,18.0,0.0,0.0),(0.0,70.0, 0.0,0.0,70.0,0.0,0.0 ,"1a")),
   "&a":((0.0,20.0,20.0,0.0, 0.0,0.0,0.0),(0.0,60.0,60.0,0.0, 0.0,0.0,0.0 ,"1a")),
   "1f":((0.0, 0.0, 0.0,0.0,47.0,0.0,0.0),(0.0, 0.0, 0.0,0.0,65.0,0.0,0.0 ,"1c")),
   "1[":((0.0, 0.0, 0.0,0.0,34.0,0.0,0.0),(0.0, 0.0, 0.0,0.0,56.0,0.0,0.0 ,"1b")),
   "4a":((0.0,40.0,40.0,0.0, 0.0,0.0,0.0),(0.0,50.0,50.0,0.0, 0.0,0.0,0.0 ,"0a")),
   "#a":((0.0,26.0,26.0,0.0, 0.0,0.0,0.0),(0.0,36.0,36.0,0.0, 0.0,0.0,0.0 ,"0a")), 
   "0i":((0.0, 0.0, 0.0,0.0,60.0,0.0,0.0),(0.0, 0.0, 0.0,0.0,60.0,0.0,0.0 ,"6n")),
   "6j":((0.0, 0.0, 0.0,0.0,60.0,0.0,0.0),(0.0, 0.0, 0.0,0.0,60.0,0.0,0.0 ,"6n"))
   }

class SuiteAngleError(Exception): pass

class SuiteAngles:
    """
    A class that stores the seven consecutive dihedral angles
    a suite consists of and allows basic calculations with them.
    """
    def __init__(self, angles):
        """
        angles is a list of seven floats.
        [deltam epsilon zeta   alpha  beta   gamma  delta]
        """
        self.set_angles(angles) # replace by property

    def __getitem__(self, i):
        return self.angles[i]

    def check_positive(self,angle):
        if angle<0: angle += 360.0
        return angle
    
    def set_angles(self, angles):
        """angles is a list of seven floats."""
        angles = map(self.check_positive,angles)
        self.deltam, self.epsilon,self.zeta,self.alpha,\
            self.beta,self.gamma,self.delta = angles
        self.angles = angles

    def confirm_angles(self):
        """checks whether all angles are in the valid range."""
        for a in self.angles:
            if not (0 <= a <= 360):
                raise SuiteAngleError('tangled: invalid angle: %6.3f'%a)

    def get_7d_diff(self, suite_angles):
        """
        Returns SuiteAngleslist with the differences
        between both suites angles.
        """
        angles = [suite_angles[i]-self[i] for i in range(7)]
        return SuiteAngles(angles)

    def dotproduct(self, suite_angles, nang):
        """
        Returns dot product of both suite angles.
        If nang=4 is set, only epsilon, zeta, alpha, beta will be used.
        Otherwise all seven
        """
        indices = nang==4 and range(1,5) or range(7)
        products = [self[i]*suite_angles[i] for i in indices]
        result = sum(products)
        while result>360: result -= 360.0
        return result

    def get_hyperellipsoid_dist(self, suite_angles, nang, weights):
        """
        Calculates distance to hyperellipsoid.
        suite contains the center point of the ellipsoid,
        weights the ellisoid radii.
        nange can be 4 or 7 to choose distance in 4D or 7D.
        Returns float
        -> Piet Hein superellipse, hyperellipsoids. 
           0      1      2      3      4      5      6  
        deltam epsilon zeta   alpha  beta   gamma  delta 
           X                                  X      X   
           X not used in 4 angle distance calc
        """
        indices = (nang==4) and range(1,5) or range(7)
        # normalize, diff < 1 inside the ellipsoid
        diffs = [ abs(self[i]-suite_angles[i])/weights[i] for i in indices ]
        powers = [ d**POWER for d in diffs ]
        return sum(powers) ** (1.0/POWER) # root

    def __str__(self):
        result = "%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f"%\
                 (self.deltam,self.epsilon,self.zeta,self.alpha,self.beta,self.gamma,self.delta)
        return result

class SuiteCluster(SuiteAngles):
    """
    Stores one of the Richardson suite clusters
    (a 7D ellipsoid).
    """
    def __init__(self,puckers, gammaname, name, \
                 dominant, satellite, wannabe, angles):
        SuiteAngles.__init__(self, angles)
        self.puckers = puckers # 33, 32, 23, 22
        self.gammaname = gammaname # p, t, or m
        self.name = name
        self.dominant = dominant # boolean
        self.satellite = satellite # boolean
        self.wannabe = wannabe # boolean


    def __str__(self):
        wannabe = self.wannabe and 'wannabe' or 'certain'
        status = 'ord'
        if self.dominant: status = 'dom'
        if self.satellite: status = 'sat'
        result = "%s\t%s%s\t%s\t%s\t"%\
                 (self.name, self.puckers,self.gammaname, status, wannabe)
        result += SuiteAngles.__str__(self)
        return result

    def get_ellipsoid_widths(self, weigh_satellites=True):
        """
        returns a 7D array of weights for the ellipsoid distance calculation.
        The weights are ellipsoid-half-widths in the given dimension.
        For: [deltam, epsilon,zeta,alpha,beta,gamma,delta]
        """
        if self.satellite and weigh_satellites:
            return SATELLITE_WEIGHTS
        else:
            return NORMAL_WEIGHTS

    def get_empirical_weights(self, domw, satw):
       """
       returns a 7D array of empirical weights for calculating
       distances between dominant and satellite clusters.
       Valid only for satellites.
       For: [deltam, epsilon,zeta,alpha,beta,gamma,delta]
       """
       # empirical weight correction for dom-sat pairs to edit
       satvalues, domvalues = DOM_SAT_EMPIRICAL[self.name]
       for i in range(7):
          satw[i] = satvalues[i] or satw[i]
          domw[i] = domvalues[i] or domw[i]
       return domw,satw


    def get_dom_sat_dist(self, dominant, suite_angles, weigh_satellites):
        """
        Recalc distances to dominant cluster and to satellite cluster
        considering different weights for points between both clusters.
        Valid only for satellite clusters.
        Returns a dist_dom, dist_sat tuple of floats.
        """
        domw = dominant.get_ellipsoid_widths(False)
        satw = self.get_ellipsoid_widths(weigh_satellites)
        domw,satw = self.get_empirical_weights(domw,satw)
        disttodom = dominant.get_hyperellipsoid_dist(suite_angles,4,False,domw)
        disttosat = self.get_hyperellipsoid_dist(suite_angles,4,False,satw)
        return disttodom, disttosat


    def get_hyperellipsoid_dist(self, suite_angles, nang, \
                                weigh_satellites=True, coordw=None):
        if not coordw:
            coordw = self.get_ellipsoid_widths(weigh_satellites)
        return SuiteAngles.get_hyperellipsoid_dist(self,suite_angles,nang,coordw)

    def is_between_dom_sat(self, dominant, suite_angles):
        """
        Called for a satellite cluster. Determines whether
        a point is in between a dominant and satellite cluster.
        if 4D dotproducts both positive, then its inbetween.
        """
        dom_suite = dominant.get_7d_diff(suite_angles)
        sat_suite = self.get_7d_diff(suite_angles)
        dom_sat = dominant.get_7d_diff(self)
        sat_dom = self.get_7d_diff(dominant)
        if (dom_suite.dotproduct(dom_sat,4) > 0) \
           and (sat_suite.dotproduct(sat_dom,4) > 0):
            return True


class ClusterSet(dict):
    """
    Specialized dictionary that reads the file with torsion angles.
    """
    def __init__(self,filename):
        '''Reads file with torsion angles.'''
        dict.__init__(self)
        for line in open(filename):
            if not line.startswith('#'):
                # parse the line
                t = line.strip().split()
                puckers, gammaname, name, unused, wannabe, color, status = t[:7]
                wannabe = wannabe=='wannabe'
                dominant = status[:3]=='dom'
                satellite = status[:3]=='sat'
                angles = [float(a) for a in t[7:]]
                # create clusters
                cluster = SuiteCluster(puckers, gammaname, name, \
                        dominant, satellite, wannabe, angles)
                # assign clusters to bins
                bin_name = puckers + gammaname
                self.setdefault(bin_name,[])
                self[bin_name].append(cluster)

    def get_bins(self):
        result = self.keys()
        result.sort()
        return result

if __name__ == '__main__':
    cl = ClusterSet('../data/suite_clusters.txt')
    for bin in cl.get_bins():
        for angles in cl[bin]:
            print angles
            #print angles.beta
        
        




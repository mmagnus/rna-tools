#!/usr/bin/python
# -*- coding: cp1252 -*-
#
# suitename.py
#
# Assigns suite names to RNA structures.
#

__author__ = "Kristian Rother, Raphael Bauer"
__credits__ = ["Marcin Domagalski","Magdalena Musielak", "Janusz Bujnicki", "Marie Curie"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Kristian Rother"
__email__ = "krother@rubor.de"
__status__ = "Production"

from suite import Suite, SuiteError, SuiteResidueError, \
     SuiteIncompleteResidueError, SuiteBfactorError,\
     SuiteDistanceError
from suite_clusters import ClusterSet
from math import pi,cos
from Scientific.Geometry import *
from dihedral import dihedral
from PDB import PDBParser
import os,sys,re

WEIGH_SATELLITES = True
WEIGH_DOM_SAT = False # True
BFACTOR_LIMIT = 60

ALLOW_WANNABE = True
CHECK_ANGLES = True
CHECK_BFACTOR = False
CHECK_DISTANCES = False

ALLOWED_ANGLES = dict(
       delta3min  =  55.0, delta3max  = 110.0,
       delta2min  = 120.0, delta2max  = 175.0,
       gammapmin  =  20.0, gammapmax  =  95.0,
       gammatmin  = 140.0, gammatmax  = 215.0,
       gammammin  = 260.0, gammammax  = 335.0,
       alphamin   =  25.0, alphamax   = 335.0,
       betamin    =  50.0, betamax    = 290.0,
       zetamin    =  25.0, zetamax    = 335.0,
       epsilonmin = 155.0, epsilonmax = 310.0
       )

"""
Automated assignment of modular suite names

(copied from the Richardson 2008 papers methods section.)

The above determination of conformer clusters was done by a
consensus of manual examination and evaluations, aided by a
variety of software for smoothing, sorting, and displaying the
multidimensional data. Now that those definitions are on hand,
there is need for an automated algorithm that can closely
approximate those suite name assignments given the list of
consensus clusters and the specific dihedral-angle values for a
new structure. The software developed to do that is called
Suitename, and the input dihedrals can be provided by Dangle
(both available from http://kinemage.biochem.duke.edu). The
dihedral-space cluster sizes and shapes are anisotropic, vary
greatly, and are mostly not Gaussian; many include only a small
number of data points; and some of the cluster pairs change
conformational roles at a boundary whose location is not readily
predictable by formula. Therefore, at the current stage of data
and understanding, the assignment algorithm and parameters are
primarily chosen to fit the manual consensus clusters rather than
determined by underlying theory.
For the general case, outer cluster boundaries are treated as
axially oriented ellipsoids, generous enough to include any
plausible data points without danger of entering an entirely
different bin. Each ellipsoid is the same size and shape, but is
centered on the mean of its cluster. The ellipsoid semi-axis in each
coordinate direction is taken as 3*sigma + 15°
, where sigma is the
average of all cluster standard deviations in that dimension. The
first term effectively scales each dimension according to its typical
range of variation, while the constant term allows for underlying
measurement uncertainties. The ellipsoid half-widths are 28° in delta,
35° in gamma, 50° in alpha, 55° in zeta, 60° in epsilon, and 70° in beta.
As noted above, data points segregate into delta and gamma bins nearly
independently of other angles, so the first step of the algorithm is
to place the data point in one of the 12 d(i-1)dg groups or else
declare it an outlier. Two ranges are accepted in d: 55°–110° (C39)
or 120°–175° (C29), and three in g: 20°–95° (p), 140°–215° (t), or
260°–335° (m), jointly defining the 12 groups. d values near zero
have been found to signify incorrect ribose stereochemistry (at
C39 for dr0004/1ET4 a212 and dr0010/1NTB b112 in the RNA05
data set), while values of e outside the range 155°–310° are usually
found to signify a misfit sugar pucker; both of these cases are
noted specifically by Suitename, as well as being named outliers.
Outliers are also declared for b outside 50°–290° or for a or z
outside 25°–335°. These single-angle outliers are referred to as
triaged.
Given the data point’s d(i-1)dg group, the next step operates in
the remaining four dimensions e, z, a, and b, to find all clusters of
which the data point could potentially be a member (that is, for
which it lies within that cluster’s ellipsoid); there may be zero, one,
or several. The scaled four-dimensional distance of the data point
from a cluster mean is zero at the mean and 1.0 anywhere on the
surface of the ellipsoid. In our initial version it was computed as a
Euclidean distance, with each component normalized by the semiaxis
in that dimension. In most cases, a data point belongs to the
cluster to which it is nearest by scaled 4D distance; this is
equivalent to drawing a boundary plane halfway between cluster
means, perpendicular to the line joining them in the scaled ezab
space. The point was given that nearest cluster’s suite name.
In order to capture more of the clearly positive manual
assignments, we found it desirable for Suitename to include
points more generously near the diagonal directions than along
the coordinate axes. This is done using a superellipsoid, the
multidimensional generalization of the superellipse (Gielis 2003)
or Lame´ oval (Gridgeman 1970). This figure bulges smoothly
outward progressively more from an ellipse into the corners of the
superscribed rectangle as the exponent increases from 2. Suitename
uses an exponent (n) of 3, in the superellipsoid equation:

|epsilon/a|^n + |zeta/b|^n + |alpha/c|^n + |beta/d|^n = 1

where a, b, c, and d are the half-widths for the relevant dihedral
angles. The data point clusters, even in the clash and B-filtered
data, show significant diagonal spread (usually in more than just
two dimensions, but not always the same ones) for highly
populated clusters such as 1a, 1c, or 1g, an effect presumably
produced either by real correlated motions or by correlated errors.
Adding the superellipsoid functionality to Suitename is still an
approximation to the probable form of the distribution in
dihedral-angle space, but it significantly improves the coverage
of what seem to be the genuine cluster boundaries.
A few very close cluster groupings require an additional
modification as well: Five dominant clusters (1a, 1c, 1b, 0a, and
6n), each in a different d(i-1)dg group, have satellite clusters with
more than half-overlapped superellipsoids. The estimated boundary
plane that was found to divide conformational types can lie
as much as four times farther from the dominant than from the
satellite cluster mean (e.g., for 1a vs. 1L). Point membership
between such pairs is decided by comparing 4D distances further
scaled in the relevant dimension (or occasionally two dimensions)
by the same ratio as of the two distances from the cluster means to
the boundary plane. If a data point is potentially a member of
more than two clusters, first its closest nondominant cluster is
found using the standard-case algorithm with default scalings, and
then the asymmetric comparison is made with the dominant
cluster if there is one. Each non-outlier suite in the input structure
is thus assigned a modular name.
To evaluate the match quality of a data point to its assigned conformational
cluster, a scaled 7D superellipsoid distance (analogous
to the 4D formula above) is now computed in all seven dihedral
dimensions of the suite, including d(i-1), g, and d as well as e, z, a,
and b. That distance d is then converted into a ‘‘suiteness’’ match
quality s, according to s = ((cos pi*d) +1) /2, which varies sinusoidally
from 1.0 at the cluster mean to zero at the surface of the
superellipsoid. A floor of 0.01 is imposed so all non-outlier points
have non-zero suiteness. All outliers have suiteness=0. The suiteness
is a measure of how well the detailed local backbone conformation
fits one of the most commonly observed (and thus presumably
most favorable) conformational clusters.
"""

class SuiteError(Exception): pass
class SuiteNoCandidateError(SuiteError): pass
class Suite7DOutlierError(SuiteError): pass
class SuiteTriageError(SuiteError): pass

class SuiteName(Suite):

   def __init__(self):
       Suite.__init__(self)
       self.message = ''

       # suite properties       
       self.puckerdm = '' # i minus 1 - pucker
       self.puckerd = '' # i - pucker
       self.gammaname = ''
       self.ddg_bin = None # string, e.g. '33p'

       self.assigned = None # assigned cluster
       self.distance = 0.0 # distance to assigned cluster
       self.suiteness=0.0

   def __str__(self):
      result = ''
      if self.assigned:
         result += "%s  assigned\n(dist %6.3f   suiteness %6.3f) \n"%(self.assigned.name,self.distance,self.suiteness)
      result += Suite.__str__(self) + '\n'
      result += self.message +'\n'
      return result
       
   def get_pucker(self, delta):
      if ALLOWED_ANGLES['delta3min'] <= delta <= ALLOWED_ANGLES['delta3max']:
         return '3'
      elif ALLOWED_ANGLES['delta2min'] <= delta <= ALLOWED_ANGLES['delta2max']:
         return '2'

   def filter_epsilon(self):
      """Checks whether epsilon is in the allowed range."""
      if self.epsilon < ALLOWED_ANGLES['epsilonmin'] or self.epsilon > ALLOWED_ANGLES['epsilonmax']:
         raise SuiteTriageError("epsilon outlier. ")

   def filter_deltam(self):
      """Checks whether delta-1 is in the allowed range."""
      self.puckerdm = self.get_pucker(self.deltam)
      if not self.puckerdm:
        raise SuiteTriageError("bad deltam. ")

   def filter_delta(self):
      """Checks whether delta is in the allowed range."""
      self.puckerd = self.get_pucker(self.delta)
      if not self.puckerd:
        raise SuiteTriageError("bad delta. ")

   def filter_gamma(self):
      """Checks whether gamma is in the allowed range."""
      if(ALLOWED_ANGLES['gammapmin'] <= self.gamma <= ALLOWED_ANGLES['gammapmax']):
           self.gammaname = 'p'
      elif(ALLOWED_ANGLES['gammatmin'] <= self.gamma <= ALLOWED_ANGLES['gammatmax']):
           self.gammaname = 't'
      elif(ALLOWED_ANGLES['gammammin'] <= self.gamma <= ALLOWED_ANGLES['gammammax']):
           self.gammaname = 'm'
      else:
           raise SuiteTriageError("gamma outlier. ")

   def filter_alpha(self):
      """Checks whether alpha is in the allowed range."""
      if(ALLOWED_ANGLES['alphamin'] > self.alpha or self.alpha > ALLOWED_ANGLES['alphamax']):
         raise SuiteTriageError("alpha outlier. ")

   def filter_beta(self):
      """Checks whether beta is in the allowed range."""
      if(ALLOWED_ANGLES['betamin'] > self.beta or self.beta > ALLOWED_ANGLES['betamax']):
         raise SuiteTriageError("beta outlier. ")

   def filter_zeta(self):
      """Checks whether zeta is in the allowed range."""
      if(ALLOWED_ANGLES['zetamin'] > self.zeta or self.zeta > ALLOWED_ANGLES['zetamax']):
         raise SuiteTriageError("zeta outlier. ")
              
   def determine_ddg_bin(self):
      """
      Returns the name of one of the 12 delta-delta-gamma bins,
      like 33p, 23t, 22m etc.
      that are used to classify suites. See Figure 2 in the
      Richardson 2008 paper.
      """
      self.filter_epsilon()
      self.filter_deltam()
      self.filter_delta()
      self.message += "pucker: %s%s. "%(self.puckerdm,self.puckerd)
      self.filter_gamma()
      self.filter_alpha()
      self.filter_beta()
      self.filter_zeta()
        
      # assign to one of 12 delta,delta,gamma bins
      self.ddg_bin = "%s%s%s"%(self.puckerdm,self.puckerd,self.gammaname)
      self.message += "ddg==%s. "%self.ddg_bin

   def get_candidate_clusters(self, cluster_set, allow_wannabe):
      """
      Returns a list of (distance, cluster) tuples for the bin
      of this suite, sorted by distance. The distance is the scaled
      normalized 4D distance of the suite to the cluster
      """
      candidates = []
      for cluster in cluster_set.get_clusters_for_bin(self.ddg_bin):
         # skip wannabe clusters
         if allow_wannabe or not cluster.wannabe:
            candidates.append((self.get_dist4d(cluster),cluster))            
      candidates.sort() # sort by distance
      for c in candidates:
         self.message += '%s(%4.2f)  '%(c[1].name,c[0])
      return candidates

   def get_closest_candidates(self, candidates):
      """
      Returns the first and second closest clusters from
      a list of (dist, cluster) tuples.
      """
      candidates += [(None,None),(None,None)]
      closest = candidates[0]
      nextclosest = candidates[1]
      if closest[1]:
         self.message += "closest: %s(%7.3f). "%(closest[1].name, closest[0])
      if nextclosest[1]:
         self.message += "nextclosest: %s(%7.3f). "%(nextclosest[1].name, nextclosest[0])
      return closest[1], nextclosest[1]

   def set_assigned_cluster(self, cluster, remark):
      """Sets this suite to the given cluster."""
      self.assigned = cluster
      self.distance = self.get_dist4d(cluster)
      self.message += remark

   def get_dominant_candidate(self,candidates):
      """Returns a dominant cluster from a list of (dist, cluster) tuples."""
      dominant = filter(lambda dc:dc.dominant, candidates)
      if dominant:
         return dominant[0]

   def get_best_satellite(self, matches):
      satellites = filter(lambda x:x.satellite, matches)
      if satellites:
         return satellites[0]

   def get_dist4d(self,cluster):
      return cluster.get_hyperellipsoid_dist(self,4, WEIGH_SATELLITES)
      
         
   def resolve_dominant_cluster(self, matches, dominant):
      """
      Match dominant cluster and closest other cluster.
      In a bin with a dominant and at least one other close candidate.
      """
      # distinguish best match cluster as satellite or just another cluster
      closest = matches[0]
      satellite = self.get_best_satellite(matches)
      if satellite:
         if satellite.is_between_dom_sat(dominant, self):
            dist_dom, dist_sat = satellite.get_dom_sat_dist(dominant,self,WEIGH_DOM_SAT)
            message = "%i-between-dom-sat (%s|%s) (%7.3f|%7.3f)"%(len(matches),dominant.name,satellite.name,dist_dom,dist_sat)
            if(dist_sat <= dist_dom):
               self.set_assigned_cluster(satellite, message)
            else:
               self.set_assigned_cluster(dominant, message)
         else:
            # point not between the means of dom and sat
            if self.get_dist4d(dominant) < self.get_dist4d(satellite):
               self.set_assigned_cluster(dominant, "%i-OUTSIDE-dominant. "%len(matches))
            else:
               self.set_assigned_cluster(satellite, "%i-OUTSIDE-satellite"%len(matches))
      else:
         # near cluster is not a satellite - take closest.
         self.set_assigned_cluster(closest,"%i-not-sat"%len(matches))

   def get_close_matches(self, candidates):
      """
      Returns a list of clusters closer than 1.0
      from a list of (distance, cluster) tuples.
      """
      cand = filter(lambda dist_c:dist_c[0]<1.0, candidates)
      return map(lambda x:x[1], cand)
   
   def determine_membership(self, cluster_set, allow_wannabe=True):
      """
      cluster membership in assigned bin find all clusters in this bin
      that match this suite conformation.
      Requires the ddg_bin to be calculated before.
      Takes into account bins that have a major cluster with satellites.

      cluster_set is a ClusterSet instance.
      allow_wannabe=True also considers wannabe clusters as candidates.
      """
      candidates = self.get_candidate_clusters(cluster_set, allow_wannabe)
      matches = self.get_close_matches(candidates)
      closest, nextclosest = self.get_closest_candidates(candidates)
      # now assign suite to a particular cluster
      if len(matches) == 1:
         # only one candidate close enough.
         self.set_assigned_cluster(closest, "1-only-one. ")      
      elif len(matches) > 1:
         # suite close to more than one cluster
         dominant = self.get_dominant_candidate(matches)
         if not dominant:
            self.set_assigned_cluster(closest,"%i-no dominant cluster."%len(matches))
         else:
            self.resolve_dominant_cluster(matches, dominant)
      elif closest:
         self.set_assigned_cluster(closest, "outlier. ")
      else:
         raise SuiteNoCandidateError('no candidates in bin %s'%self.ddg_bin)
      self.calc_suiteness(closest)

   def calc_suiteness(self,closest):
      """recompute distance for all 7 dimensions, calc suiteness"""
      self.suiteness = 0.0
      if self.assigned:
         dist = self.assigned.get_hyperellipsoid_dist(self,7, WEIGH_SATELLITES)
         if dist > 1.0:
            message = "7D distance forces assignment of %s(%6.3f) to be outlier"%(self.assigned.name,dist)
            message += '\n'+str(self)+'\n'
            message += str(self.assigned)
            self.assigned = None
            raise Suite7DOutlierError(message)
      else:
         dist = closest.get_hyperellipsoid_dist(self,7,False)
         if dist <= 1.0:
            self.set_assigned_cluster((closest,dist),"7D distance forces assignment to closest cluster: %s(%6.3f)"%(self.assigned.name,dist))

      if (dist <= 1): # not an outlier
         self.suiteness = (cos(pi*dist) +1)/2.0
         self.suiteness = max(self.suiteness, 0.01) # floor
      
   def assign_name(self, cs, allow_wannabe=ALLOW_WANNABE, check_angles=CHECK_ANGLES,\
      check_bfactor=CHECK_BFACTOR, check_distances=CHECK_DISTANCES):
      """
      Returns the name of the suite, if it can be assigned.
      """
      try:
         self.confirm(check_angles, check_bfactor, check_distances)
         self.determine_ddg_bin()
         self.determine_membership(cs, allow_wannabe)
         if self.assigned:
            return self.assigned.name
      except SuiteError,e:
         self.message += str(e)
         raise e

class PDBProcessor:
   def __init__(self):
      self.good = 0
      self.resierrors = 0
      
      self.counts = {}
      self.out = []

   def get_string_for_chain(self,chain):
      pass
      
   def get_strings_for_pdbfile(self,pdb_file_name):
      parser = PDBParser()
      structure= parser.get_structure(pdb_file_name.replace('.pdb',''),pdb_file_name)
      for chain in structure[0].child_list:
         chain_id=chain.id
         result = self.get_string_for_chain(chain)
         if not re.search('^-+$',result):
            yield chain_id,result

   def count(self,seq):
      for i in range(0,len(seq),2):
         code = seq[i:i+2]
         self.counts.setdefault(code,0)
         self.counts[code] += 1


   def loop_path(self,path):
      for fn in os.listdir(path):
         if fn[-4:] in ['.pdb','.ent','.PDB','.ENT']:
            for chain_id,suite_string in self.get_strings_for_pdbfile(path + fn):
               outline = path+fn+'\t0\t'+chain_id+'\t'+suite_string
               print outline
               self.out.append(outline+'\n')
               self.count(suite_string)

   def report(self):
      pass

   def write_output(self,outfile):
      open(outfile,'w').writelines(self.out)


class SuiteNameProcessor(PDBProcessor):

   def __init__(self):
      PDBProcessor.__init__(self)
      self.cs = ClusterSet()
      self.triaged = 0
      self.outlier7d = 0
      self.emptybin = 0

   def report(self):
      keys = self.counts.keys()
      keys.sort()
      for k in keys:
         print k,self.counts[k]
      print "suites assigned",snp.good
      print "bad residues   ",snp.resierrors
      print "triaged resi   ",snp.triaged
      print "7D outliers    ",snp.outlier7d
      print "empty bins     ",snp.emptybin

   def get_string_for_chain(self,chain):
      result = ''
      first = None
      second = None
      for resi in chain.child_list:
         first, second = second, resi
         if first and second:
            try:
               suite = SuiteName()
               suite.set_residues(first, second)
               name = suite.assign_name(self.cs)
               self.good += 1
               result += name
            except SuiteTriageError, e:
               self.triaged += 1
               result += 'tt'
            except Suite7DOutlierError, e:
               self.outlier7d += 1
               result += 'oo'
            except SuiteNoCandidateError, e:
               self.emptybin += 1
               result += '!!'
            except SuiteResidueError, e:
               self.resierrors += 1
               result += '--'
      return result
               

if __name__ == '__main__':
   print """
erena - Suite Suffix Tree
(c) 2008 Kristian Rother and Raphael Bauer
Python implementation of the suitename program by Jane Richardson et al.

usage: suitename.py <path with pdb files> <outfile>
   """   
   snp = SuiteNameProcessor()
   path = sys.argv[1]
   outfile = sys.argv[2]
   snp.loop_path(path)
   snp.report()
   snp.write_output(outfile)

   
   


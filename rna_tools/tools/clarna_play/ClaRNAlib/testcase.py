import os
import re
import unittest
import numpy as np

class TestClassifier(unittest.TestCase):

    def setUp(self):
        self.PYTHON_BIN="/opt/local/bin/python"
        self.CL_LIB="/tmp/classifier-2012-07-31.json.gz"
        
    def test_doublets_class(self):
        from utils import PDBObject,DoubletsDict,GraphTool,bench_start,bench_stop
        from distances import doublet_params_dict, Doublet, Residue, residue_conformation

        eps = 0.001

        for pdb_id in ['3fo6','2zjp']:
            dd = DoubletsDict(reduced_atoms=['*'])
            dd.load_pdb(pdb_id)
            gr = GraphTool(PDBObject.pdb_fn(pdb_id,"close_doublets"),edge_type="dist")
            
            old_params = {}
            bench_start("old params computations")
            for d_id in gr.get_ids():
                full_d_id = pdb_id+":"+d_id
                n_type = dd.get_n_type(full_d_id)
                (p1,p2) = dd.get(full_d_id)
                for type in ['bp','stacking','base-phosphate','base-ribose']:
                    params = doublet_params_dict((p1,p2), n_type, type)
                    old_params["%s-%s"%(d_id,type)] = params
            bench_stop("old params computations")
            
            residues = {}
            bench_start("new params computations")
            for d_id in gr.get_ids():
                full_d_id = pdb_id+":"+d_id
                n_type = dd.get_n_type(full_d_id)
                # print "processing doublet: %s (%s)" % (full_d_id,n_type)
                (p1,p2) = dd.get(full_d_id)
                (r1,r2) = d_id.split(":")
                if not residues.has_key(r1):
                    residues[r1] = Residue(r1,n_type[0],p1)
                if not residues.has_key(r2):
                    residues[r2] = Residue(r2,n_type[1],p2)
                d = Doublet(d_id,residues[r1],residues[r2])
                for type in ['bp','stacking','base-phosphate','base-ribose']:
                    params = old_params["%s-%s"%(d_id,type)]
                    for key,expected_value in params.items():
                        if re.match("^(dist_[A-Z]|i_|ii_|oxygens)",key):
                            if type=='base-phosphate':
                                v = d.ph_info.get(key)
                            elif type=='base-ribose':
                                v = d.br_info.get(key)
                            else:
                                v = None
                        else:
                            v = getattr(d, key)
                        # print "d_id=%(d_id)s type=%(type)s key=%(key)s v=%(v)s expected=%(expected_value)s" % locals()
                        if key in ['oxygens']:
                            self.assertTrue(sorted(v)==sorted(expected_value))
                        else:
                            self.assertTrue(v>expected_value-eps and v<expected_value+eps,"%s: got: %.4f, expected: %.4f"%(key,v,expected_value))
            bench_stop("new params computations")

    def test_distances(self):
        # 1A9N:R5:R17 -- WW_cis
        # 1M90:A379:A386 -- WW_cis
        # 2UXC:A394:A476 -- WH_cis
        from utils import DoubletsDict
        from distances import bp_distance, bp_distance2
        r_atoms = ['N1','N2','N3','N4','N6','N7','P']
        r_atoms += ['C2','C4','C5','C6','C8',"C1'"]
        r_atoms += ['O2']
        r_atoms += ["O2'","O3'","O4'"]
        r_atoms += ["OP1","OP2","O5'","NEXT:O3'"]
        r_atoms += ['N1','C6','O6','C5','C4','N3','C2','N2','N7','C8','N9']

        dd = DoubletsDict(reduced_atoms=r_atoms)
        doublets = ['1A9N:R5:R17','1M90:A379:A386','2UXC:A394:A476']
        print "distance between %s and %s: %.5f" % (doublets[0], doublets[1], bp_distance(dd.get(doublets[0]), dd.get(doublets[1]), dd.get_n_type(doublets[0])))
        print "distance between %s and %s: %.5f" % (doublets[0], doublets[2], bp_distance(dd.get(doublets[0]), dd.get(doublets[2]), dd.get_n_type(doublets[0])))
        print "distance2 between %s and %s: %.5f" % (doublets[0], doublets[1], bp_distance2(dd.get(doublets[0]), dd.get(doublets[1]), dd.get_n_type(doublets[0])))
        print "distance2 between %s and %s: %.5f" % (doublets[0], doublets[2], bp_distance2(dd.get(doublets[0]), dd.get(doublets[2]), dd.get_n_type(doublets[0])))
    
    def test_bp_params(self):
        from utils import DoubletsDict
        from distances import doublet_params_dict, Doublet, Residue
        
        TEST_DATA = [
            (
                "3FO6:A39:A57",
                {
                    'stack_orient': 3, 
                    'dist_z': 0.22861798604329428, 
                    'o_ang': 69.303115470181339, 
                    'dist': 5.5976623073175231, 
                    'stack_overlap': 0, 
                    'conf': 0, 
                    'n12cc_ang': 87.910385850994587, 
                    'n1cc_ang': 92.089614149005413, 
                    'strand_orient_norm': -1, 
                    'stack_norm': 0.97095424, 
                    'nn_ang_norm': 13.698485786637832, 
                    'nn_ang': 166.30151421336217, 
                    'min_dist': 1.7859186, 
                    'strand_orient': -1, 
                    'orient': 1, 
                    'n2cc_ang': 88.640047368908554
                }
            ),
            (
                "2ZJP:X303:X77",
                {
                    'stack_orient': 3, 
                    'strand_orient': -1, 
                    'dist': 6.3820344572213887, 
                    'stack_min_dist': 1.4906554379247356, 
                    'n12cc_ang': 79.337562855790168, 
                    'stack_norm': 0.93983656, 
                    'strand_orient_norm': -1, 
                    'conf': 0, 
                    'n2_z': -0.93983656, 
                    'orient': 1, 
                    'n2cc_ang': 96.415019366401566, 
                    'dist_z': 1.1562830607096353, 
                    'o_ang': 50.822836557417922, 
                    'stack_overlap': 0, 
                    'rot_ang': -56.801295513103156, 
                    'n1cc_ang': 100.66243714420983, 
                    'min_dist': 1.4906554379247356, 
                    'nn_ang_norm': 20.020604505050017, 
                    'nn_ang': 159.97939549494998
                }
            )
        ]
        eps = 0.001
        
        dd = DoubletsDict(reduced_atoms=['*'])
        
        for d_id, expected_params in TEST_DATA:
            n_type = dd.get_n_type(d_id)
            points = dd.get(d_id)
            self.assertEqual(len(points),2)
            params = doublet_params_dict(points, n_type, 'stacking')
            self.assertTrue(params is not None)
    
            for key,expected_value in expected_params.items():
                self.assertTrue(params.has_key(key))
                self.assertTrue(params[key]>expected_value-eps) 
                self.assertTrue(params[key]<expected_value+eps) 
            
            # test nowych metod
            r1 = Residue("A1",n_type[0],points[0])
            r2 = Residue("B1",n_type[1],points[1])
            d = Doublet(d_id,r1,r2)
            for key,expected_value in expected_params.items():
                if key in ['dist','min_dist','nn_ang','nn_ang_norm','n1cc_ang','n2cc_ang',
                    'n12cc_ang','o_ang','orient','stack_orient','stack_norm','strand_orient',
                    'strand_orient_norm','conf','dist_z','n2_z','rot_ang','stack_min_dist',
                    'stack_overlap']:
                    v = getattr(d,key)
                    # print "%s got: %.4f, expected: %.4f" % (key,v,expected_value)
                    self.assertTrue(v>expected_value-eps and v<expected_value+eps,"%s: got: %.4f, expected: %.4f"%(key,v,expected_value))
                else:
                    print "skipping: %s" % key


# workflows
automatix xx or makex

```shell

[mm] evox$ git:(master) ✗ xxdownload_job_and_process.py ade ade_cp15-ca11f0d6
mkdir: simrna_ade_cp15-ca11f0d6: File exists
/usr/local/lib/python2.7/site-packages/urllib3/connectionpool.py:858: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
  InsecureRequestWarning)
['ade_cp15', 'ca11f0d6', 'thrs6.50A_clust04X.
```

```shell

set -x
rna_pdb_tools.py --un-nmr *_MD.pdb
for i in `ls *MD_*`; do rna_pdb_tools.py --rpr $i --replace-hetatm  --dont-fix-missing-atoms > ${i}_rpr.pdb; done	
for i in `ls *rpr*`; do rna_pdb_tools.py --edit 'A:1-35>A:51-85,B:1-11>B:20-30' $i > ${i}_fixChains.pdb; done
```

```shell
cat test.sh
# X
# rna_pdb_tools.py --un-nmr yC_5LJ3_U2U6_core_mdrFx_addh_MD.pdb
set -x
#for i in `ls *MD_*` do rna_pdb_tools.py --rpr $i --replace-hetatm  --dont-fix-missing-atoms > ${i}_rpr.pdb done

for i in `ls *rpr*`; do
	rna_pdb_tools.py --edit 'A:1-35>A:51-85,B:1-11>B:20-30' $i > ${i}_fixChains.pdb
done

```

```shell

for i in `ls -d */`;
	 do
	     echo $i;
	     cd $i
	     rm *.pdb
	     cp ../1st-yC_5LJ3_U2U6_core_mdr_chainAB+Extra2ResiAtEnds_delt1-1+dC85_orient.pdb .
	     cp ../test.sh .
	     ./test.sh & # to make it semi-parallel ;-)
	     cd ..
	 done;
```

Minimal workflow system.

See for more robust workflow management system:

- https://snakemake.github.io
- https://www.nextflow.io

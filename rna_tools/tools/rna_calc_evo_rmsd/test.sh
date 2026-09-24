#./rna_calc_evo_rmsd.py -a test_data/aln.sto -t test_data/4qk8_cl.pdb test_data/4qlm_cl.pdb -n 4qk8_cl -m test_data/mapping.txt
./rna_calc_evo_rmsd.py -a test_data/aln.sto -t test_data/4qk8_cl.pdb test_data/4qlm_cl.pdb

./rna_calc_evo_rmsd.py -a test_data/aln_nox.sto -t test_data/4qk8_cl.pdb test_data/4qlm_cl.pdb -v

./rna_calc_evo_rmsd.py -a test_data/trna.sto -t test_data/1ehz_std.pdb test_data/6Y2L_2_std.pdb

# automated, no alignment needed: Rfam (requires Infernal and Rfam.cm, see RfamAlign.py)
#./rna_calc_evo_rmsd.py --rfam_db Rfam.cm -t test_data/1ehz_std.pdb test_data/6Y2L_2_std.pdb -v
# automated, pairwise sequence alignment
./rna_calc_evo_rmsd.py --auto pairwise -t test_data/1ehz_std.pdb test_data/6Y2L_2_std.pdb
./rna_calc_evo_rmsd.py --auto pairwise -t test_data/4qk8_cl.pdb test_data/4qlm_cl.pdb -v

# Rfam family example, RF01750 (ZMP/ZTP riboswitch): get the list of 3D structures of the family
# from Rfam, download them (one PDB file per chain), map them on Rfam, align them to the family
# covariance model and calculate RMSD and sequence identity (esl-alipid, Easel) for all pairs;
# plot: test_data/RF01750/RF01750_rmsd_seqid_vs_rmsd.png. Requires internet, Infernal, Easel and biopython.
FAM=RF01750
D=test_data/$FAM
mkdir -p $D
# structures of the family: rfam_acc pdb_id chain pdb_start pdb_end bit_score evalue_score cm_start cm_end hex_colour
curl -sS https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.pdb.gz | gunzip -c | awk -v fam=$FAM '$1 == fam' > $D/${FAM}_structures.tsv
echo "# structures of $FAM in Rfam:"
cat $D/${FAM}_structures.tsv
# covariance model of the family, used as the Rfam db (instead of the whole Rfam.cm)
curl -sS -o $D/$FAM.cm https://rfam.org/family/$FAM/cm
cmpress -F $D/$FAM.cm > /dev/null
# download structures (mmCIF, works also for big entries) and save each chain as <pdb>_<chain>.pdb
awk '{print tolower($2), $3}' $D/${FAM}_structures.tsv | sort -u | while read pdb chain; do
    [ -f $D/${pdb}.cif ] || curl -sS -o $D/${pdb}.cif https://files.rcsb.org/download/${pdb}.cif
    python -c "
import sys
from Bio.PDB import MMCIFParser, PDBIO, Select
cif, chain, out = sys.argv[1:]
s = MMCIFParser(QUIET=True).get_structure('', cif)
class ChainSelect(Select):
    def accept_model(self, m): return m.id == 0
    def accept_chain(self, c): return c.id == chain
io = PDBIO(); io.set_structure(s); io.save(out, ChainSelect())
" $D/${pdb}.cif $chain $D/${pdb}_${chain}.pdb
done
# all vs all, alignment saved to test_data/RF01750/RF01750_rmsd_rfam/RF01750.sto
./rna_calc_evo_rmsd.py --all_vs_all --rfam_db $D/$FAM.cm -o $D/${FAM}_rmsd.csv $D/*_*.pdb

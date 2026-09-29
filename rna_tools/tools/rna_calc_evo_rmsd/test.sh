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
    def accept_residue(self, r): return r.id[0] == ' '  # skip waters/ions/ligands (some use
        # 5-char CCD codes that don't fit the legacy PDB format and break parsers)
io = PDBIO(); io.set_structure(s); io.save(out, ChainSelect())
" $D/${pdb}.cif $chain $D/${pdb}_${chain}.pdb
done
# align every extracted chain to the family CM (independent of the --all_vs_all cut_ga
# filtering below) just to show the structures table with sequences aligned to the model
python -c "
from RfamAlign import write_fasta, run
from RNAmodel import RNAmodel
from rna_calc_evo_rmsd import get_sequence
from RNAalignment import RNAalignment

rows = []
seqs = []
with open('$D/${FAM}_structures.tsv') as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        pdb, chain = cols[1].lower(), cols[2]
        name = '%s_%s' % (pdb, chain)
        seq = get_sequence(RNAmodel.parse('$D/%s.pdb' % name))
        rows.append((cols, name))
        seqs.append((name, seq))

fasta = '$D/${FAM}_structures.fa'
aln_fn = '$D/${FAM}_structures.sto'
write_fasta(seqs, fasta)
run(['cmalign', '-o', aln_fn, '$D/$FAM.cm', fasta])
ra = RNAalignment(aln_fn)

print('# structures of $FAM in Rfam, aligned to the CM:')
for cols, name in rows:
    print('\t'.join(cols + [ra.get_seq(name)]))
"
# all vs all, alignment saved to test_data/RF01750/RF01750_rmsd_rfam/RF01750.sto
# --rfam_evalue 0.01 matches Infernal's default significance ('!') threshold instead of
# the curated Rfam gathering threshold (--cut_ga), so borderline hits like 5btp_A are kept
./rna_calc_evo_rmsd.py --all_vs_all --rfam_db $D/$FAM.cm --rfam_evalue 0.01 -o $D/${FAM}_rmsd.csv $D/*_*.pdb

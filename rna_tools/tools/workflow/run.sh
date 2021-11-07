for i in ` find . -iname 'struc'`;
do
    cd $i
    parallel /Users/magnus/work/src/rna-tools/rna_tools/tools/mq/rsRNASP/run_rsRNASP ::: *.pdb # > rsRNASP
    cd /Volumes/Toshiba/mq/decoys/_others_/rasp/rasp
done

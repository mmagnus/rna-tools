folder=job`openssl rand -hex 12`
ssh mq mkdir $folder
rsync --delete --include '*.pdb' --exclude="*" -rv $1/* mq:$folder
# rna_mq_collect.py
ssh mq "cd $folder && /home/mqapRNA/mqaprna_env/mqapRNA/main/mqaprna * -o mq.csv"
rsync -r  --exclude='*.pdb'  --exclude='.csvsort*'  mq:$folder/ $1

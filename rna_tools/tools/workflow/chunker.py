


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

#main
f = open('pdb_list.txt').read().split()
# 30 = runs * 10 = 300cores
no_of_files = int(round(len(f) / 100))
c = 0
for i in chunks(f, no_of_files):
    txt = "echo '" 
    # -g farna_hires.csv
    txt += "/home/magnus/mqaprna_env/mqapRNA/main/mqaprna  -o FARNA__hires_" + str(c) + " " + " ".join(i) +"'"
    txt += " | qsub -N mq" + str(c) + " -V -cwd -pe mpi 1"
    print txt
    c += 1
    

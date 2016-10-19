## RUN:
## Rscript bpseq2cl.R csv_file.csv
## input: csv_file.csv
## output: csv_file.csv_CLnotation.txt

# These 3 lines are to call on console the file for input
args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
bpseq <- read.csv(fn, header = FALSE) 
names(bpseq) <- c("posi", "res", "posj")
bpseq

# Creating the new data frame object 
df <- data.frame(chainA=character(nrow(bpseq)),
                 posi=integer(nrow(bpseq)),
                 chain=character(nrow(bpseq)),
                 posj=integer(nrow(bpseq)),
                 bptype=character(nrow(bpseq)),
                 residuei=character(nrow(bpseq)),
                 residuej=character(nrow(bpseq)),
                 typeligation=character(nrow(bpseq)),
                 value=numeric(nrow(bpseq)))
df$chainA <- rep('A',nrow(bpseq))
df$chain <- rep('A',nrow(bpseq))
df$bptype <- rep('bp',nrow(bpseq))
df$typeligation <- rep('WW_cis',nrow(bpseq))
df$value <- rep(round(1, digits=2),nrow(bpseq))
df$posi <- bpseq$posi
df$posj <- bpseq$posj
df$residuei <- bpseq$res
df$residuei <- as.character(df$residuei)

# To keep only the resi and resj that actually form a connection
df <- df[!df$posj == 0, ] 
df$residuej <- bpseq$res[bpseq$posj[1:nrow(bpseq)]]
df$residuej <- as.character(df$residuej)

# To remove repeated connections. Exemple: resi=6 and resj=378 are connected, as are resi=378 with resj=6. We don't need the repetition.
cols = c(2,4) 
newdf = df[,cols]
for (i in 1:nrow(df))
{
  newdf[i, ] = sort(df[i,cols])
}
df <- df[!duplicated(newdf),]

# clarna_compare.py has as first line "Classifier: Clarna", so this line must be appended
df <- rbind(c("Classifier: Clarna", "", "", "", "", "", "", "", ""), df)
df

write.table(df, file = paste(args[1],"_CLnotation.txt", sep = ""), col.names = F, row.names = F, sep = " ")

# Before running clarna_compare.py on your files, remove all "
args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
bpseq <- read.csv(fn, header = FALSE)
names(bpseq) <- c("posi", "res", "posj")
bpseq

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

df <- df[!df$posj == 0, ]

df$residuej <- bpseq$res[bpseq$posj[1:nrow(bpseq)]]
#df
df <- df[1:(nrow(df)/2),]
df

write.table(df, file = paste(args[1],"_CLnotation.txt", sep = ""), col.names = F, row.names = F, sep = "\t")

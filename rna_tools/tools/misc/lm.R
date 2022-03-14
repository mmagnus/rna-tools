f = "/Users/magnus/work/src/rna-tools/rna_tools/tools/triplexibility/t12.csv"
df <- read.table(f, header=TRUE, sep=';')
df
head(df)
df
l = lm(growth ~ ., data=df)
#l = glm(growth ~ ., data=df, family = binomial)
summary(l)
#plot(df$t1_3_rmsd ~ df$growth)
 
cor(df)

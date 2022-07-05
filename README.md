# MATS: Example

Here, we provide a quick example using the MATS package for a continuous outcome.

### Simulate Genotypes
First, we simulate genotype data from 2 populations (EUR and AFR) using the sim1000G package.

```
library("sim1000G")

ped_file_1000genomes = system.file("examples", "20130606_g1k.ped", package = "sim1000G")
ped = read.table(ped_file_1000genomes,h=T,as=T,sep="\t")
examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region-chr4-93-TMEM156.vcf.gz" )

pop1 = c("TSI", "FIN", "GBR", "IBS")
id1 = ped$Individual.ID [ ped$Population %in% pop1 ]

pop2 = c("YRI", "LWK", "MAG", "MSL", "ENA", "ASW", "ACB")
id2 = ped$Individual.ID [ ped$Population %in% pop2 ]


vcf  = readVCF(vcf_file, min_maf = 0.02, max_maf = 0.32)

startSimulation(vcf, subset = c(id1))
n1 = 300
id.EUR = c()
for(i in 1:n1) id.EUR[i] = SIM$addUnrelatedIndividual()
genotypes.EUR = retrieveGenotypes(id.EUR)
rownames(genotypes.EUR) = paste0("EUR", 1:n1)

startSimulation(vcf, subset = c(id2))
n2 = 300
id.AFR = c()
for(i in 1:n2) id.AFR[i] = SIM$addUnrelatedIndividual()
genotypes.AFR = retrieveGenotypes(id.AFR)
rownames(genotypes.AFR) = paste0("AFR", 1:n2)

p = 10
genotypes = rbind(genotypes.EUR, genotypes.AFR)[,1:p]

```

### Set Remaining Parameters

```
n = nrow(genotypes)
Y = rnorm(n)
trait_type = "continuous"

w = matrix(rnorm(p), ncol = 1)
xhat = genotypes %*% w
groups = c(rep("EUR", n1), rep("AFR", n2))

Sex=sample(c("F", "M"), n, replace = TRUE)
Age=as.numeric(floor(runif(n, 30, 100)))
C = as.matrix(data.frame(Sex, Age))
categorical.vars = c("Sex")

Y = 0.05 + 0.1*xhat + xhat*ifelse(groups == "AFR", 0, 0.1) + rnorm(n) 

nP = 2
P <- prcomp(genotypes, center = TRUE, scale = TRUE)$x
ev <- prcomp(genotypes, center = TRUE, scale = TRUE)$sdev^2
K <- cov(t(scale(genotypes, center = TRUE, scale = TRUE)))

pc_main_effects = TRUE
```

# Run MATS

```
MATS(Y, xhat, groups, C, P, ev, K, categorical.vars, trait_type, np, pc_main_effects)
```





# A guide to running the code
## Description of each file and its role in the order they are run

1. [**getAllSNPs.R**](/getAllSNPs.R) - gets a list of random SNP IDs available from SNPedia
2. [**getSNPNames.py**](/getSNPNames.py) - Uses a list of completely randomly chosen SNP IDs and picks a set number of IDs per chromosome to output to a list of the curated SNPs
3. [**getRiskAlleles.py**](/getRiskAlleles.py) - From a list of SNPs known to predispose CD, keeps the SNPs with sufficient data in Ensembl
4. [**frequencyEthnicities.py**](/frequencyEthnicities.py) - For every single SNP chosen, gathers its allele frequency data from 1000GP and pools by ancestral region
5. [**conductRegression.py**](/conductRegression.py) - For each SNP, conducts linear regression analysis against two separate variables and ouputs results to new file
6. [**smallerRSquareCD.py**](/smallerRSquareCD.py) - Looking through the regression results for SNPs against CD prevalence. For each risk SNP, calculates the percentage of non-risk SNPs for which it has a greater association
7. [**smallerRSquareWheat.py**](/smallerRSquareWheat.py) - Looking through the regression results for SNPs against wheat agriculture history. For each risk SNP, calculates the percentage of non-risk SNPs for which it has a greater association

## Dependencies
- [statsmodels](https://www.statsmodels.org/stable/index.html#)

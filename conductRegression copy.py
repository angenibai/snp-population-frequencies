import statsmodels.api as sm
import csv

# The CD Prevalences and Wheat History Categories for each ethnicity in order
cdPrevalence = [0.67,0.96,0,0.03,1.29,0.64,0.76,0]
wheatHist = [1,1,4,3,2,3,2,4]

# Header for regression summaries
regressionHead = ["SNP","coefficient","R-squared","Adj R-squared","p-value"]

with open("pooled_snps copy.csv", "r") as readF:
    read = csv.reader(readF)
    dataList = list(read)

# Setting up to write to two different files
writeWheat = [regressionHead]
writeCD = [regressionHead]

for row in dataList[1:]:
    snp = row[0]
    alleleFreq = list(map(float, row[4:]))

    wheatHistConst = sm.add_constant(wheatHist)
    alleleFreqConst = sm.add_constant(alleleFreq)

    # Independent wheat history and dependent allele frequency
    result = sm.OLS(alleleFreq, wheatHistConst).fit()
    writeWheat.append([snp,result.params[1],result.rsquared,result.rsquared_adj,result.pvalues[1]])

    # Independent allele frequency and dependent CD prevalence
    result = sm.OLS(cdPrevalence, alleleFreqConst).fit()
    writeCD.append([snp,result.params[1],result.rsquared,result.rsquared_adj,result.pvalues[1]])

with open("regression_results_wheat copy.csv","w") as writeF:
    write = csv.writer(writeF)
    write.writerows(writeWheat)

with open("regression_results_cd copy.csv","w") as writeF:
    write = csv.writer(writeF)
    write.writerows(writeCD)

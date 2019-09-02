import statsmodels.api as sm
import csv

"""
For each SNP, conducts linear regression analysis against two separate variables
and ouputs results to new file
"""

"""
SET FILENAMES BELOW
"""
popDataFile = "snp_population_data.csv"
resultsCDFile = "regression_results_wheat.csv"
resultsWheatFile = "regression_results_cd.csv"

"""
START INDEX
May start from different points in popDataFile but remember that the results
are output to a new file - which WILL rewrite previous files with same name
"""
startIdx = 695

"""
Data for CD Prevalences and Wheat Agriculture Categories for each ethinicty
in order
"""
cdPrevalence = [0.67,1.05,0,0.03,1.29,0.64,0.76,0]
wheatHist = [1,1,4,3,2,3,2,4]

# Header for regression summaries
regressionHead = ["SNP","Risk","coefficient","R-squared","Adj R-squared","p-value"]
popsBegin = 4

with open(popDataFile, "r") as readF:
    read = csv.reader(readF)
    dataList = list(read)

# Setting up to write to two different files
writeWheat = [regressionHead]
writeCD = [regressionHead]

for row in dataList[startIdx:]:
    snp = row[0]
    risk = row[3]
    alleleFreq = list(map(float, row[popsBegin:]))

    wheatHistConst = sm.add_constant(wheatHist)
    alleleFreqConst = sm.add_constant(alleleFreq)

    # Independent wheat history and dependent allele frequency
    result = sm.OLS(alleleFreq, wheatHistConst).fit()
    writeWheat.append([snp,risk,result.params[1],result.rsquared,result.rsquared_adj,result.pvalues[1]])

    # Independent allele frequency and dependent CD prevalence
    result = sm.OLS(cdPrevalence, alleleFreqConst).fit()
    writeCD.append([snp,risk,result.params[1],result.rsquared,result.rsquared_adj,result.pvalues[1]])

with open(resultsWheatFile,"w") as writeF:
    write = csv.writer(writeF)
    write.writerows(writeWheat)

with open(resultsCDFile,"w") as writeF:
    write = csv.writer(writeF)
    write.writerows(writeCD)

import csv

"""
Looking through the regression results for SNPs against wheat agriculture history.
For each risk SNP, calculates the percentage of non-risk SNPs for which it
has a greater association
"""

"""
SET FILENAMES BELOW
"""
regressionResultFile = "regression_results_wheat.csv"
calculatedResultFile = "smaller_rsquares_wheat.csv"

rSquareIdx = 3
pValIdx = 5

with open(regressionResultFile, "r") as readF:
    read = csv.reader(readF)
    readList = list(read)
    nonRisk = [row for row in readList if row[1] == "None"]
    risk = [row for row in readList if row[1] == "Non-HLA" or row[1] == "HLA"]

nonRisk.sort(key=lambda x:x[rSquareIdx])

for riskRow in risk:
    riskSNP = riskRow[0]
    riskType = riskRow[1]
    riskR = float(riskRow[rSquareIdx])

    lessThan = 0
    for nonRiskRow in nonRisk:
        curR = float(nonRiskRow[rSquareIdx])
        if curR < riskR:
            lessThan+=1
        else:
            break

    with open(calculatedResultFile, "r") as readF:
        read = csv.reader(readF)
        dataList = list(read)

    dataList.append([riskSNP, riskType, riskR, lessThan, float(lessThan/652*100), riskRow[pValIdx]])

    with open(calculatedResultFile, "w") as writeF:
        csvWriter = csv.writer(writeF)
        csvWriter.writerows(dataList)

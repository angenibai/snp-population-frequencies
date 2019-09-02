import sys, csv, requests

"""
From a list of SNPs known to predispose CD, keeps the SNPs with sufficient
data in Ensembl
"""

"""
SET FILENAMES
"""
readAllRiskFile = "CD risk alleles.csv"
chosenRiskFile = "CD risk alleles.csv"

with open(readAllRiskFile, "r") as readF:
    read = csv.reader(readF)
    dataList = list(read)

dataWrite = [["SNP","Trait","Risk Allele"]]

for row in dataList[1:]:
    snp = row[0]
    print(snp)
    # Grabbing data about phenotypes
    server = "https://rest.ensembl.org"
    ext = "/variation/human/"+snp+"?phenotypes=1"
    r = requests.get(server+ext, headers={"Content-Type" : "application/json"})

    if not r.ok:
        continue

    rawData = r.json()

    phenotypes = rawData["phenotypes"]

    # Making sure we actually have phenotypes
    if not phenotypes:
        continue

    # Finding the celiac disease one
    for p in phenotypes:
        if p["trait"] == "CELIAC DISEASE":
            if "risk_allele" in p.keys():
                risk = p["risk_allele"]
                # Checking that the risk allele is actually mapped
                mapping = rawData["mappings"][0]
                if risk in mapping["allele_string"]:
                    newRow = [snp, p["trait"], p["risk_allele"]]
                    print("Risk allele found ", p["risk_allele"])
                    dataWrite.append(newRow)
            break

with open(chosenRiskFile, "w") as writeF:
    csvWriter = csv.writer(writeF)
    csvWriter.writerows(dataWrite)

import requests, sys, csv

targetAlleles = {
    "rs3184504":"T",
    "rs2327832":"G",
    "rs13151961":"A"
}

'''
European --> 0
Mediterranean --> 1
Native American --> 2
Northeast Asian --> 3
Northern European --> 4
Southeast Asian --> 5
Southwest Asian --> 6
Sub-Saharan African --> 7
'''

EURO = 0
MEDI = 1
AMER = 2
NEAN = 3
NEUR = 4
SEAN = 5
SWAN = 6
SAFR = 7
OTHR = 8
popLen = 8

popToEth = { "CHB":NEAN, "JPT":NEAN, "CHS":NEAN, "CDX":NEAN, "KHV":NEAN,
    "CEU":NEUR, "TSI":MEDI, "FIN":NEUR, "GBR":NEUR, "IBS":MEDI, "YRI":SAFR,
    "LWK":SAFR, "GWD":SAFR, "MSL":SAFR, "ESN":SAFR, "ASW":SAFR, "ACB":OTHR,
    "MXL":OTHR, "PUR":MEDI, "CLM":EURO, "PEL":AMER, "GIH":SWAN, "PJL":SWAN,
    "BEB":SEAN, "STU":SEAN, "ITU":SWAN
}

ethNames = ["European", "Mediterranean", "Native American", "Northeast Asian", "Northern European",
"Southeast Asian", "Southwest Asian", "Sub-Saharan African"]

unwantedPops = ["AFR","AMR","EAS","EUR","SAS","ACB","MXL","ALL"]

for snp in targetAlleles:
    target = targetAlleles[snp]
    print("Collecting data for %s target %s" %(snp, target))

    # Calling data part
    server = "https://rest.ensembl.org"
    ext = "/variation/human/"+snp+"?pops=1"

    r = requests.get(server+ext, headers={"Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    rawData = r.json()

    # Initialising allele counts and frequencies for this snp
    alleleCounts = [0 for i in range(popLen)]
    alleleFreqs = [0.0 for i in range(popLen)]
    totalPop = [0 for i in range(popLen)]

    for p in rawData["populations"]:
        name = p["population"]
        if "1000GENOMES:phase_3" not in name:
            # We only care about populations from 1000 Genomes Phase 3
            continue
        pop = name[len(name)-3:]
        if pop in unwantedPops:
            # Again filtering out the ones we don't want
            continue

        # Grabbing info
        eth = popToEth[pop]
        newCount = p["allele_count"]
        newFreq = p["frequency"]

        if p["allele"] == target:
            # We can update our data for this population's ethnicity

            totalPop[eth] += int(newCount/newFreq + 0.5)
            alleleFreqs[eth] = (alleleCounts[eth]+newCount)/totalPop[eth]

            alleleCounts[eth] += newCount

        elif p["frequency"] == 1:
            # When we know that there is zero of the target allele
            totalPop[eth] += int(newCount/newFreq + 0.5)
            alleleFreqs[eth] = alleleCounts[eth]/(totalPop[eth])

    # Time to put this into the csv file
    with open("pooled_snps.csv","r",encoding="utf-8") as readF:
        read = csv.reader(readF)
        dataList = list(read)

    dataList[0].append("%s (%s)" %(snp, target)) # Adding the name of the SNP as a header
    for eth in range(popLen):
        # Adding frequency of this SNP onto the end of the row for this ethnicity
        dataList[eth+1].append(alleleFreqs[eth])

    with open("pooled_snps.csv","w",encoding="utf-8") as writeF:
        csvWriter = csv.writer(writeF)
        csvWriter.writerows(dataList)

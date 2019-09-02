import requests, sys, csv
from time import sleep

"""
For every single SNP chosen, gathers its allele frequency data from 1000GP
and pools by ancestral region
"""

"""
SET FILENAMES
"""
readSNPFile = "chosen_snps.csv"
popDataFile = "snp_population_data.csv"

"""
SET START INDEX
This program adds onto existing files
"""
startIdx = 708

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

with open(readSNPFile, "r") as readF:
    read = csv.reader(readF)
    snpList = list(read)

counter = 1
for row in snpList[708:]:
    snp, target, chr, risk = row
    print("%d Collecting data for %s target %s" %(counter, snp, target))
    counter += 1

    # Gotta make sure we don't go over the request limit
    if counter % 50 == 0:
        sleep(45)

    # Calling data part
    server = "https://rest.ensembl.org"
    ext = "/variation/human/"+snp+"?pops=1"

    r = requests.get(server+ext, headers={"Content-Type" : "application/json"})

    if not r.ok:
        continue

    rawData = r.json()

    '''
    # Check minor allele frequency is acceptable
    if not rawData["MAF"]:
        continue
    if rawData["MAF"] < 0.01:
        continue
    '''

    # Get some basic info
    m = rawData["mappings"][0]
    chr = m["seq_region_name"]
    location = m["location"]

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
            alleleFreqs[eth] = float(alleleCounts[eth]+newCount)/totalPop[eth]
            alleleCounts[eth] += newCount

        elif p["frequency"] == 1:
            # When we know that there is zero of the target allele
            totalPop[eth] += int(newCount/newFreq + 0.5)
            alleleFreqs[eth] = alleleCounts[eth]/(totalPop[eth])


    # Adding to an existing population file
    with open(popDataFile,"r") as readF:
        read = csv.reader(readF)
        dataList = list(read)

    # Name, chr, location, target allele, then ethnicities in order
    newRow = [snp,target,chr,risk] + alleleFreqs
    dataList.append(newRow)

    with open(popDataFile,"w") as writeF:
        csvWriter = csv.writer(writeF)
        csvWriter.writerows(dataList)

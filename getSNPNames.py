import requests, sys, csv

"""
Uses a list of completely randomly chosen SNP IDs and picks a set number of
IDs per chromosome to output to a list of the curated SNPs
"""

"""
SET FILENAMES
"""
readSNPFile = "snpList-8000.csv"
chosenSNPFile = "chosen_snps.csv"

"""
SET MAXIMUM NUMBER OF SNPS PER CHROMOSOME
"""
MAX_PER_CHR = 30

"""
SET MINIMUM BASE PAIRS OF SEPARATION BETWEEN SNPS
"""
minSeparation = 3000

# The three bits of data we are interested in for each snp
snps = [[] for i in range(23)]
majorAlleles = [[] for i in range(23)]
occupiedPositions = [[] for i in range(23)]

# Checks to make sure we don't have too many
totalFull = 0
full = [False for i in range(23)]

with open(readSNPFile,"r") as readF:
    read = csv.reader(readF)
    dataList = list(read)

for row in dataList[1:]:
    snp, c, p = row #  c and p are pretty unreliable data
    snp = snp.lower()

    # Server calling for information
    server = "https://rest.ensembl.org"
    ext = "/variation/human/"+snp+"?"
    r = requests.get(server+ext, headers={"Content-Type" : "application/json"})

    if not r.ok:
        continue

    rawData = r.json()

    # We want the mapping with a numerical chromosome
    print(snp)
    chrom = "X"
    for m in rawData["mappings"]:
        try:
            # Getting some basic info
            chrom = int(m["seq_region_name"])
            pos = m["start"]
            break
        except ValueError:
            continue

    # Making sure we have a numerical chromosomes
    try:
        int(chrom)
    except ValueError:
        #print("Rejected because chromosome isn't valid")
        continue
    if chrom > 22 or chrom < 1:
        #print("Rejected because chromosome isn't numerical: ", chrom)
        continue

    # Checking if this chromosome is already full
    if full[chrom]:
        #print("Chromsome %d is already full" %(chrom))
        continue

    # Check minor allele frequency is acceptable
    maf = rawData["MAF"]
    if not maf:
        #print("No minor allele frequency")
        continue
    if maf < 0.01:
        #print("MAF < 0.01: ", maf)
        continue

    majorAllele = rawData["ancestral_allele"]
    if not majorAllele:
        #print("No major allele")
        continue

    # Check it's not too close to anything
    tooClose = False
    occupied = occupiedPositions[chrom]
    if occupied:
        for p in occupied:
            if abs(p - pos) < minSeparation:
                tooClose = True
                break
    if tooClose:
        #print("It's too close to something")
        continue

    # Ok it's passed all our tests
    snps[chrom].append(snp)
    majorAlleles[chrom].append(majorAllele)
    occupiedPositions[chrom].append(pos)

    print("Chromsome %d filled SNP number %d" %(chrom, len(snps[chrom])))

    # Check if this chromosome is filled with snps
    if len(snps[chrom]) == MAX_PER_CHR:
        full[chrom] = True
        totalFull += 1

    # Check if we're totally done
    if totalFull == 22:
        break

if totalFull < 22:
    print("We are not full yet")

# Time to write to our own csv file!
dataWrite = [["SNP","Target Allele","Chromosome"]]
for c in range(1,23):
    for s in range(len(snps[c])):
        newRow = [snps[c][s], majorAlleles[c][s], c]
        dataWrite.append(newRow)

with open(chosenSNPFile, "w") as writeF:
    csvWriter = csv.writer(writeF)
    csvWriter.writerows(dataWrite)

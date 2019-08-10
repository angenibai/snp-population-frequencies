import requests, sys, csv

# The three bits of data we are interested in for each snp
snps = [[] for i in range(23)]
majorAlleles = [[] for i in range(23)]
occupiedPositions = [[] for i in range(23)]

# Checks to make sure we don't have too many
MAX_PER_CHR = 10
totalFull = 0
full = [False for i in range(23)]

with open("snpList copy.csv","r") as readF:
    read = csv.reader(readF)
    dataList = list(read)

for row in dataList[1:]:
    snp, chr, pos = row #  chr and pos are pretty unreliable data
    snp = snp.lower()

    # Server calling for information
    server = "https://rest.ensembl.org"
    ext = "/variation/human/"+snp+"?"
    r = requests.get(server+ext, headers={"Content-Type" : "application/json"})

    if not r.ok:
        continue

    rawData = r.json()

    # Getting some basic info
    m = rawData["mappings"][0]
    chr = m["seq_region_name"]
    pos = m["start"]

    # Making sure we have a numerical chromosomes
    try:
        chr = int(chr)
    except ValueError:
        continue
    if chr > 22 or chr < 1:
        continue

    # Checking if this chromosome is already full
    if full[chr]:
        continue

    # Check minor allele frequency is acceptable
    maf = rawData["MAF"]
    if not maf:
        continue
    if maf < 0.01:
        continue

    majorAllele = rawData["ancestral_allele"]
    if not majorAllele:
        continue

    # Check it's not too close to anything
    tooClose = False
    occupied = occupiedPositions[chr]
    if occupied:
        for p in occupied:
            if abs(p - pos) < 3000:
                tooClose = True
                break
    if tooClose:
        continue

    # Ok it's passed all our tests
    snps[chr].append(snp)
    majorAlleles[chr].append(majorAllele)
    occupiedPositions[chr].append(pos)

    print("Chromsome %d filled SNP number %d" %(chr, len(snps[chr])))

    # Check if this chromosome is filled with snps
    if len(snps[chr]) == MAX_PER_CHR:
        full[chr] = True
        totalFull += 1

    # Check if we're totally done
    if totalFull == 22:
        break


if totalFull < 22:
    print("We are not full yet")

# Time to write to our own csv file!
dataWrite = [["SNP","Target Allele","Chromosome"]]
for c in range(1,23):
    for s in range(MAX_PER_CHR):
        newRow = [snps[c][s], majorAlleles[c][s], c]
        dataWrite.append(newRow)

with open("chosen_snps copy.csv", "w") as writeF:
    csvWriter = csv.writer(writeF)
    csvWriter.writerows(dataWrite)

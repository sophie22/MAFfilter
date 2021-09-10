#Python 3
import sys

import pandas as pd
import vcf

from CellBase_api import build_url, get_response

# Store the clinically relevant transcripts in a Python list object
# Obtain list of clinically relevant transcripts from the genes2transcript (g2t) text file
g2t_list = []
with open("MAFfilter/g2tr.txt", 'r') as f:
    for line in f:
        g2t_list.append(line.strip().split('\t')[1])
# print(g2t_list[:10])

# Initialise the dataframe to store information and annotation of the filtered variants
variants = pd.DataFrame(columns=["chromosome", "position", "ref", "alt", "qual", "MAF", "HGNC", "RefSeq", "HGVS", "consequence"])

def transcript_check(csq_list):
    """
        check if any of the transcripts annotated for the variant
        are in the genes2transcript list, if so return the annotation
        if not, return the annotations for every NM_ transcript
    """
    # csq_transcripts = [i.split('|')[3] for i in csq_list]
    found = False
    for i in csq_list:
        csq = i.split('|')
        if csq[3] in g2t_list:
            found = True
            return [csq[2],csq[3],csq[12],csq[5]] # "HGNC", "RefSeq", "HGVS", "consequence"
    if found is False:
        csqs = [i.split('|') for i in csq_list]
        NMs = [nm for nm in csqs if nm[3].startswith("NM_")]
        # "HGNC", "RefSeq", "HGVS", "consequence"
        HGNCs = ",".join([h[2] for h in NMs])
        RefSeqs = ",".join([n[3] for n in NMs])
        HGVSs = ",".join([g[12] for g in NMs])
        cons = ",".join([c[5] for c in NMs])
        return [HGNCs, RefSeqs,HGVSs,cons]

# Open the patient vcf that was filtered by bcftools to only contain homozygous variants
inp = sys.argv[1]
vcf_reader = vcf.Reader(open(inp, 'r'))

# Read in the information about the homozygous variants
for index, record in enumerate(vcf_reader):
    # print(index)
# for record in vcf_reader:
    chrom = str(record.CHROM)
    pos = str(record.POS)
    ref = record.REF
    alt = str(record.ALT[0])
    qual = str(record.QUAL)
    # create the variant identifier to be used to query the CellBase API
    v = "%3A".join([chrom,pos,ref,alt])
    url = build_url(v)
    data_result = get_response(url)
    # parse the MAF from the CellBase response
    try:
        popFreqStudies = data_result["populationFrequencies"]
        try: # obtain MAF from the gnomAD exomes all ethnicities field (gnEall)
            gnEall = [i["altAlleleFreq"] for i in popFreqStudies if i["study"] == "GNOMAD_EXOMES" and i["population"] == "ALL"]
            if gnEall[0] < 0.01:
                csq_list = record.INFO['CSQ'] # check whether transcript is in the g2t list and populate annotation fields
                check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
                # add the variant information and annotation to the dataframe as it has MAF < 0.01
                variants = variants.append({
                    "chromosome": chrom, "position": pos, "ref": ref, "alt": alt, "qual": qual,
                    "MAF": gnEall[0], "HGNC": check[0], "RefSeq": check[1], "HGVS": check[2],
                    "consequence": check[3]
                }, ignore_index=True)
        except: # obtain MAF from the gnomAD genomes all ethnicities field (gnGall)
            gnGall = [j["altAlleleFreq"] for j in popFreqStudies if j["study"] == "GNOMAD_GENOMES" and j["population"] == "ALL"]
            if gnGall[0] < 0.01:
                csq_list = record.INFO['CSQ'] # check whether transcript is in the g2t list and populate annotation fields
                check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
                # add the variant information and annotation to the dataframe as it has MAF < 0.01
                variants = variants.append({
                    "chromosome": chrom, "position": pos, "ref": ref, "alt": alt, "qual": qual,
                    "MAF": gnGall[0], "HGNC": check[0], "RefSeq": check[1], "HGVS": check[2],
                    "consequence": check[3]
                }, ignore_index=True)
    # if there is no gnomAD exomes or genomes allele frequency information available,
    # assume that variant is absent from gnomAD and has MAF of 0
    except:
        csq_list = record.INFO['CSQ'] # check whether transcript is in the g2t list and populate annotation fields
        check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
        # add the variant information and annotation to the dataframe as it has MAF < 0.01
        variants = variants.append({
            "chromosome": chrom, "position": pos, "ref": ref, "alt": alt, "qual": qual,
            "MAF": 0, "HGNC": check[0], "RefSeq": check[1], "HGVS": check[2],
            "consequence": check[3]
        }, ignore_index=True)

# Write variant information and annotation from the dataframe to a tsv file
variants.to_csv("X111923_hom_MAFeg0.01_annot.tsv", sep='\t', encoding='utf-8', index=False)
print("Done!")

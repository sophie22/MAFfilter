#vcf parser 3
import sys
import vcf
import pandas as pd
from CallBase_api import build_url, get_response

g2t_list = []
with open("MAFfilter/g2tr.txt", 'r') as f:
    for line in f:
        g2t_list.append(line.strip().split('\t')[1])
# print(g2t_list[:10])

inp = sys.argv[1]
variants = pd.DataFrame(columns=["chromosome", "position", "ref", "alt", "MAF", "HGVS", "HGNC", "RefSeq", "consequence", "other"])
vcf_reader = vcf.Reader(open(inp, 'r'))

def transcript_check(csq_list):
    """
    check if any of the transcripts annotated for the variant
    are in the GEMINI genes2transcript list
    """
    # csq_transcripts = [i.split('|')[3] for i in csq_list]
    found = False
    for i in csq_list:
        csq = i.split('|')
        if csq[3] in g2t_list:
            found = True
            return [csq[1],csq[2],csq[3],csq[5],[""]] # "HGVS", "HGNC", "RefSeq", "consequence", "other"
    if found is False:
        csq_transcripts = [j.split('|')[3] for j in csq_list]
        NM_transcripts = [k for k in csq_transcripts if k.startswith("NM_")]
        return ["","","","",NM_transcripts]

for record in vcf_reader:
    # record = next(vcf_reader)
    chrom = str(record.CHROM)
    pos = str(record.POS)
    ref = record.REF
    alt = str(record.ALT[0])
    v = ":".join([chrom,pos,ref,alt])
    url = build_url(v)
    data_result = get_response(url)
    try:
        popFreqStudies = data_result["populationFrequencies"]
        gnEall = [i["altAlleleFreq"] for i in popFreqStudies if i["study"] == "GNOMAD_EXOMES" and i["population"] == "ALL"]
        if gnEall[0] < 1:
            csq_list = record.INFO['CSQ']
            check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
            variants = variants.append({
                "chromosome": chrom, "position": pos, "ref": ref, "alt": alt,
                "MAF": gnEall[0], "HGVS": check[0], "HGNC": check[1],
                "RefSeq": check[2], "consequence": check[3], "other": check[4]
            }, ignore_index=True)

    except:
        pass
        # print("") #data not available

# print(variants.head(15))
variants.to_csv("X111923_annot_MAF_filtered.tsv", sep='\t', encoding='utf-8', index=False)
print("Done!")
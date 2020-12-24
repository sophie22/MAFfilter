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
variants = pd.DataFrame(columns=["chromosome", "position", "ref", "alt", "qual", "MAF", "HGNC", "RefSeq", "HGVS", "consequence"])
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
            # "HGNC", "RefSeq", "HGVS", "consequence"
            return [csq[2],csq[3],csq[12],csq[5]]
    if found is False:
        csqs = [i.split('|') for i in csq_list]
        NMs = [nm for nm in csqs if nm[3].startswith("NM_")]
        # "HGNC", "RefSeq", "HGVS", "consequence"
        HGNCs = ",".join([h[2] for h in NMs])
        RefSeqs = ",".join([n[3] for n in NMs])
        HGVSs = ",".join([g[12] for g in NMs])
        cons = ",".join([c[5] for c in NMs])
        return [HGNCs, RefSeqs,HGVSs,cons]

for index, record in enumerate(vcf_reader):
    print(index)
# for record in vcf_reader:
    chrom = str(record.CHROM)
    pos = str(record.POS)
    ref = record.REF
    alt = str(record.ALT[0])
    qual = str(record.QUAL)
    v = "%3A".join([chrom,pos,ref,alt])
    url = build_url(v)
    data_result = get_response(url)
    try:
        popFreqStudies = data_result["populationFrequencies"]
        try:
            gnEall = [i["altAlleleFreq"] for i in popFreqStudies if i["study"] == "GNOMAD_EXOMES" and i["population"] == "ALL"]
            if gnEall[0] < 0.01:
                csq_list = record.INFO['CSQ']
                check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
                variants = variants.append({
                    "chromosome": chrom, "position": pos, "ref": ref, "alt": alt, "qual": qual,
                    "MAF": gnEall[0], "HGNC": check[0], "RefSeq": check[1], "HGVS": check[2],
                    "consequence": check[3]
                }, ignore_index=True)
        except:
            gnGall = [j["altAlleleFreq"] for j in popFreqStudies if j["study"] == "GNOMAD_GENOMES" and j["population"] == "ALL"]
            if gnGall[0] < 0.01:
                csq_list = record.INFO['CSQ']
                check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
                variants = variants.append({
                    "chromosome": chrom, "position": pos, "ref": ref, "alt": alt, "qual": qual,
                    "MAF": gnGall[0], "HGNC": check[0], "RefSeq": check[1], "HGVS": check[2],
                    "consequence": check[3]
                }, ignore_index=True)

    except:
        csq_list = record.INFO['CSQ']
        check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
        variants = variants.append({
            "chromosome": chrom, "position": pos, "ref": ref, "alt": alt, "qual": qual,
            "MAF": 0, "HGNC": check[0], "RefSeq": check[1], "HGVS": check[2],
            "consequence": check[3]
        }, ignore_index=True)
        # break

variants.to_csv("X111923_hom_MAFeg0.01_annot.tsv", sep='\t', encoding='utf-8', index=False)
print("Done!")
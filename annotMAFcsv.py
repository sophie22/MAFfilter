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
variants = pd.DataFrame(columns=["chromosome", "position", "ref", "alt", "MAF", "HGNC", "RefSeq", "HGVS", "consequence", "other"])
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
            return [csq[2],csq[3],csq[12],csq[5],""] # "HGNC", "RefSeq", "HGVS", "consequence", "other"
    if found is False:
        csq_transcripts = [j.split('|')[3] for j in csq_list]
        NM = [k for k in csq_transcripts if k.startswith("NM_")]
        NM_transcripts = ",".join(NM)
        return ["","","","",NM_transcripts]

failed = []
with open("X111923_failedList2.txt", "w") as vf:
    for index, record in enumerate(vcf_reader):
        # while index < 2000:
        print(index)
        # for record in vcf_reader:
        chrom = str(record.CHROM)
        pos = str(record.POS)
        ref = record.REF
        alt = str(record.ALT[0])
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
                        "chromosome": chrom, "position": pos, "ref": ref, "alt": alt,
                        "MAF": gnEall[0], "HGNC": check[0],
                        "RefSeq": check[1], "HGVS": check[2], "consequence": check[3], "other": check[4]
                    }, ignore_index=True)
            except:
                gnGall = [j["altAlleleFreq"] for j in popFreqStudies if j["study"] == "GNOMAD_GENOMES" and j["population"] == "ALL"]
                if gnGall[0] < 0.01:
                    csq_list = record.INFO['CSQ']
                    check = transcript_check(csq_list) # "HGVS", "HGNC", "RefSeq", "consequence", "other"
                    variants = variants.append({
                        "chromosome": chrom, "position": pos, "ref": ref, "alt": alt,
                        "MAF": gnGall[0], "HGNC": check[0],
                        "RefSeq": check[1], "HGVS": check[2], "consequence": check[3], "other": check[4]
                    }, ignore_index=True)

        except:
            variant = "-".join([chrom,pos,ref,alt])
            vf.write(str(index) + "\t" + variant+"\n")
            failed.append(index)

            # break

# print(variants.head(15))
variants.to_csv("X111923_annot_MAFeg0.01_filtered2.tsv", sep='\t', encoding='utf-8', index=False)
# print(failed)
print("Done!")
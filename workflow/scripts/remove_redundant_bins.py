import os
import sys
import shutil

def redundant_bins(seqs_dir, output_dir):
    print(seqs_dir)
    vrhyme_dir = os.path.dirname(seqs_dir)
    bins_file = ""
    for filename in os.listdir(vrhyme_dir):
        if filename.endswith(".summary.tsv"):
            bins_file = (os.path.join(vrhyme_dir, filename))

    idx = 0
    list_bins_ok = []
    with open(bins_file, 'r') as f:
        for line in f:
            if idx == 0:
                idx += 1
            else:
                elmt = line.split('\t')
                if int(elmt[3]) < 2:
                    list_bins_ok.append("vRhyme_bin_"+elmt[0]+".fasta")
                    list_bins_ok.append("vRhyme_bin_"+elmt[0]+".faa")
                    list_bins_ok.append("vRhyme_bin_"+elmt[0]+".ffn")

    if not(os.path.exists(output_dir)):
        os.mkdir(output_dir)
    for item in list_bins_ok:
        shutil.copyfile(os.path.join(seqs_dir, item), os.path.join(output_dir, item))

redundant_bins(sys.argv[1], sys.argv[2])

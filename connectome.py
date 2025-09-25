#!/usr/bin/env python
from umi_tools import network
import pysam
import sys
import json
import pandas as pd
import os


### Extract 30-mer, CB and UMI from bam file and append to DataFrame
def bam2df (fin):
    samfile = pysam.AlignmentFile(fin, "rb", check_sq=False)
    Mer = []
    CB = []
    UMI = []
    for line in samfile:
        line=line.to_string()
        if line.strip():
            if len(line.split("\t"))==23:
                mer = line.split("\t")[9]
                cb = line.split("\t")[19]
                umi = line.split("\t")[22]
                #test = [line.split("\t")[9], line.split("\t")[19], line.split("\t")[22]]
                #read = "\t".join([line.split("\t")[i] for i in [9,19,22]])
                Mer.append(mer)
                CB.append(cb)
                UMI.append(umi)
    #print('CB type=', type(CB))
    df = pd.DataFrame({'Mer':Mer, 'CB':CB, 'UMI':UMI})
    samfile.close()
    #print(df.head())
    df['Mer']=df.Mer.str[20:50]
    #print(df.head())
    df=df.groupby(['Mer','CB','UMI']).size().reset_index(name='count').sort_values(by='count', ascending=False)
    return df

## First collapse of the extracted 30mer-CB-UMI to identical lines + trim to 30-mer
### Cut synbar to 30-mer and collapse tables to one occurence of "30mer-CB-UMI"
#def first_collapse (df):
#	df['Mer']=df.Mer.str[20:50]
#    df=df.groupby(['Mer','CB','UMI']).size().reset_index(name='count').sort_values(by='count', ascending=False)
#    return df


### Demuxing extracted 30-mers on a 'per cell barcode' basis  
def demux_nmer(df,fileout, hamming):
	grouped = df.groupby("CB")
	for name, group in grouped:
		cb = []
		cbc = dict()
		for row_index, row in group.iterrows():
			cols=row['Mer']+'|'+row['CB']+'|'+row['UMI']
			cols = bytes(cols, "ascii")
			cb.append(cols)
			cbc[cols]=int(row['count'])
		uc = network.UMIClusterer()
		CBclusters = uc(cbc, int(hamming))
		cbFinal = dict()
		for l in CBclusters:  # This is a list with the first element the dominant barcode
			cbFinal[l[0]] = 0
			for x in l:  # Iterate over all barcodes in a cluster
				cbFinal[l[0]] += cbc[x]
		
        ## write to a final barcode table file
		with open(fileout, 'a') as f:
			for k, v in cbFinal.items():
				k=str(k)
				k = k.replace("b", "").replace("'", "")
				k="\t".join(k.split("|"))
				f.write(k+'\t'+str(v) + '\n')
### Cut synbar to 30-mer and collapse tables to one occurence of "30mer-CB-UMI"


def main (fin, fout, hamming):
    if os.path.isfile(sys.argv[2]):
        os.remove(sys.argv[2])
    fin = sys.argv[1]
    fout = sys.argv[2]
    hamming = sys.argv[3]
    print('Extracting barcodes from BAM')
    df = bam2df(fin)
    print('Demuxing 30-mers')
    demux_nmer(df,fout, hamming)

def usage():
    print("Usage:   python ", sys.argv[0]," file_in file_out hamming_distance")
    print("'file_in' is a bam file")
    print("'file_out' is final table of collapsed barcodes") 
    print("'hamming_distance' for 30-mer dedup, usually 1")
    print('collapses identical 30mer-CB-UMI, using user-defined hamming distance')



if __name__=='__main__':
    if len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        usage()


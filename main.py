# import libraries
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import time


def add_ref_alt(file, out_file):
    df = pd.read_csv(file, encoding='latin1')

    try:
        mutation_cds = df.iloc[:, 20]
        for i in range(len(df)):
            if mutation_cds[i][-3:] == 'del':
                df.loc[df.index[i], 'ref'] = 'del'
                df.loc[df.index[i], 'alt'] = 'del'
            elif mutation_cds[i][-2].isupper():
                df.loc[df.index[i], 'ref'] = mutation_cds[i][-2]
                df.loc[df.index[i], 'alt'] = mutation_cds[i][-1]
            else:
                df.loc[df.index[i], 'ref'] = mutation_cds[i][-3]
                df.loc[df.index[i], 'alt'] = mutation_cds[i][-1]
    except IndexError:
        for i in range(len(df)):
            df.loc[df.index[i], 'ref'] = '.'
            df.loc[df.index[i], 'alt'] = '.'
        print('index error!')
    df.to_csv(out_file, index=False)


# gnomAD file
def tsv_to_csv(file, file2, name):
    # reading given tsv file
    with open(file, 'r') as myfile:
        with open(file2, 'w') as csv_file:
            for line in myfile:
                # Replace every tab with comma
                fileContent = re.sub("\t", ",", line)

                # Writing into csv file
                csv_file.write(fileContent)

    df = pd.read_csv(file2)
    df.rename(columns={'AF': name}, inplace=True)
    df.to_csv(file2)


# gnomAD file
def add_genomic_coordinate(afr, est, amr, asj, fin, eas, nwe, seu):
    files = [afr, est, amr, asj, fin, eas, nwe, seu]

    for file in files:

        df = pd.read_csv(file, dtype={'CHR': 'string', 'BP': 'int64', 'ref': 'string', 'alt': 'string', 'AC': 'int64',
                                      'AF_amr': 'float', 'AN': 'float', 'homozygote_count': 'float', 'L2': 'float'})

        chrom = df.iloc[:, 1]
        #print(chrom)
        pos = df.iloc[:, 2]
        #print(pos)
        for i in range(len(df)):
            new_pos = int(pos[i]) - 1
            if chrom[i] == 'X' or chrom[i] == 'Y':
                df.loc[df.index[i], 'MUTATION_GENOME_POSITION'] = '23:' + str(new_pos) + '-' + str(
                    new_pos)
            else:
                df.loc[df.index[i], 'MUTATION_GENOME_POSITION'] = str(chrom[i]) + ':' + str(new_pos) + '-' + str(
                    new_pos)
        df.to_csv(file, index=False)


# Merge cosmic and gnomAD
def merge_cosmic_gnomad(cosmic, file_2, file_3):
    print(f"merging {file_2}...", flush=True)

    # Pyarrow: https://pythonspeed.com/articles/pandas-read-csv-fast/
    gnomad = pd.read_csv(file_2, engine='pyarrow')

    output1 = pd.merge(cosmic, gnomad,
                       on=['MUTATION_GENOME_POSITION', 'ref', 'alt'],
                       how='inner')

    output1.to_csv(file_3)
    print(f"{file_2} completed.", flush=True)


def concat_gnomad(afr, est, amr, asj, eas, fin, nwe, seu, merge_csv):
    afr = pd.read_csv(afr)
    #print(afr)
    est = pd.read_csv(est)
    #print(est)
    amr = pd.read_csv(amr)
    #print(amr)
    asj = pd.read_csv(asj)
    #print(asj)
    eas = pd.read_csv(eas)
    #print(eas)
    fin = pd.read_csv(fin)
    #print(fin)
    nwe = pd.read_csv(nwe)
    #print(nwe)
    seu = pd.read_csv(seu)
    #print(seu)

    concat = pd.concat([afr, est, amr, asj, eas, fin, nwe, seu])

    concat.drop_duplicates(subset=['MUTATION_GENOME_POSITION','ref', 'alt', 'AF_afr', 'AF_est', 'AF_amr', 'AF_asj', 'AF_eas', 'AF_fin','AF_nwe', 'AF_seu'], inplace=True)

    concat.to_csv(merge_csv)


def make_histogram(merged_file, gene):
    # open merged cosmic and gnomAD file
    file = pd.read_csv(merged_file)

    if len(file) > 0:

        positions = list(file['MUTATION_GENOME_POSITION'])

        afr = list(file['AF_afr'])

        est = list(file['AF_est'])

        amr = list(file['AF_amr'])

        asj = list(file['AF_asj'])

        eas = list(file['AF_eas'])

        fin = list(file['AF_fin'])

        nwe = list(file['AF_nwe'])

        seu = list(file['AF_seu'])

        # Functions to make multi-bar histogram
        X = positions

        X_axis = np.arange(len(X))

        plt.figure(figsize=(40, 25))
        plt.bar(range(len(afr)), afr, align='edge', width=0.5, color='#00c5ff', label='African')
        plt.bar(range(len(est)), est, align='edge', width=0.5, color='#ff00f1', label='Estonian')
        plt.bar(range(len(amr)), amr, align='edge', width=0.5, color='#0400ff', label='American Admixed/Latino')
        plt.bar(range(len(asj)), asj, align='edge', width=0.5, color='#8a00ff', label='Ashkenazi Jewish')
        plt.bar(range(len(eas)), eas, align='edge', width=0.5, color='#ff0039', label='East Asian')
        plt.bar(range(len(fin)), fin, align='edge', width=0.5, color='#ff6500', label='Finnish European')
        plt.bar(range(len(nwe)), nwe, align='edge', width=0.5, color='#19e634', label='Non-Western European')
        plt.bar(range(len(seu)), seu, align='edge', width=0.5, color='#ffe000', label='Southern European')

        plt.xticks(X_axis, X, rotation=45, size=15)
        plt.xlabel('Positions', fontsize=25)
        plt.ylabel('Allele Frequency', fontsize=25)
        plt.title(gene, fontsize=40)
        plt.legend(loc='best', borderpad=2, prop={"size": 16})

        # save histogram
        plt.savefig(f'./results/{gene}.png', dpi=500)

        plt.close()

        # Clear lists
        positions.clear()

        afr.clear()

        est.clear()

        amr.clear()

        asj.clear()

        eas.clear()

        fin.clear()

        nwe.clear()

        seu.clear()

def main(c_file, c_file2, g_afr_tsv, g_afr_csv, g_est_tsv, g_est_csv, g_amr_tsv, g_amr_csv, g_asj_tsv, g_asj_csv,
         g_eas_tsv, g_eas_csv, g_fin_tsv, g_fin_csv, g_nwe_tsv, g_nwe_csv, g_seu_tsv, g_seu_csv, merged, gene, n_cores):
    # Cosmic
    add_ref_alt(c_file, c_file2)

    # Gnomad
    tsv_to_csv(g_afr_tsv, g_afr_csv, 'AF_afr')
    tsv_to_csv(g_est_tsv, g_est_csv, 'AF_est')
    tsv_to_csv(g_amr_tsv, g_amr_csv, 'AF_amr')
    tsv_to_csv(g_asj_tsv, g_asj_csv, 'AF_asj')
    tsv_to_csv(g_eas_tsv, g_eas_csv, 'AF_eas')
    tsv_to_csv(g_fin_tsv, g_fin_csv, 'AF_fin')
    tsv_to_csv(g_nwe_tsv, g_nwe_csv, 'AF_nwe')
    tsv_to_csv(g_seu_tsv, g_seu_csv, 'AF_seu')

    add_genomic_coordinate(g_afr_csv, g_est_csv, g_amr_csv, g_asj_csv, g_eas_csv, g_fin_csv, g_nwe_csv, g_seu_csv)

    # Merge data
    cosmic = pd.read_csv(c_file2)
    with multiprocessing.Pool(n_cores) as p:
        p.starmap(merge_cosmic_gnomad, [(cosmic, g_afr_csv, 'afr.csv'),
                                        (cosmic, g_est_csv, 'est.csv'),
                                        (cosmic, g_amr_csv, 'amr.csv'),
                                        (cosmic, g_asj_csv, 'asj.csv'),
                                        (cosmic, g_eas_csv, 'eas.csv'),
                                        (cosmic, g_fin_csv, 'fin.csv'),
                                        (cosmic, g_nwe_csv, 'nwe.csv'),
                                        (cosmic, g_seu_csv, 'seu.csv')])

    # concat all the merged files together
    concat_gnomad('afr.csv', 'est.csv', 'amr.csv', 'asj.csv', 'eas.csv', 'fin.csv', 'nwe.csv', 'seu.csv',
                  merged)

    # Make Histograms
    make_histogram(merged, gene)


#######
if __name__ == '__main__':
    N_CORES = 4

    files = os.listdir('./0')

    for i in range(len(files)):
        gene = files[i][:-4]

        start_time = time.time()
        main(f'./0/{files[i]}', 'cosmic_final.csv',
             'gnomad.genomes.r2.1.1.afr.adj.ld_scores.ldscore', 'gnomad_afr.csv',
             'gnomad.genomes.r2.1.1.est.adj.ld_scores.ldscore', 'gnomad_est.csv',
             'gnomad.genomes.r2.1.1.amr.adj.ld_scores.ldscore', 'gnomad_amr.csv',
             'gnomad.genomes.r2.1.1.asj.adj.ld_scores.ldscore', 'gnomad_asj.csv',
             'gnomad.genomes.r2.1.1.eas.adj.ld_scores.ldscore', 'gnomad_eas.csv',
             'gnomad.genomes.r2.1.1.fin.adj.ld_scores.ldscore', 'gnomad_fin.csv',
             'gnomad.genomes.r2.1.1.nwe.adj.ld_scores.ldscore', 'gnomad_nwe.csv',
             'gnomad.genomes.r2.1.1.seu.adj.ld_scores.ldscore', 'gnomad_seu.csv',
             f'./results/{gene}_merged.csv', gene, N_CORES)


        print(f"Total time: {time.time() - start_time:.2f} secs.")

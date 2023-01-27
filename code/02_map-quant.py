# Imports
import os
import glob
import logging
from collections import defaultdict

def main():
    # Init
    logging.basicConfig(level=logging.INFO,
                    filename=snakemake.log[0],
                    format='%(message)s')
                    
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info('#'*80)
    logging.info('#####                                                                      #####')
    logging.info('#####                      Mapping and quantification                      #####')
    logging.info('#####                                                      *** / *****     #####')
    logging.info('#'*80)
    logging.info('#####\n##')
    logging.info('##\tMapping and quantifying 10x libraries:')

    indexes_10x = index_list(snakemake.config['pin_10x'])
    for index, seq_names in indexes_10x.items():
        process_index(index, seq_names, '10x')

    if (snakemake.config['workflow'] == '10x + HTO'):
        indexes_hto = index_list(snakemake.config['pin_hto'])
        for index, seq_names in indexes_hto.items():
            process_index(index, seq_names, 'hto')

    os.system(f'touch {snakemake.output}')

def index_list(pin_dict):
    indexes = defaultdict(list)
    for seq_name in pin_dict:
        for index in pin_dict[seq_name]:
            indexes[index] += [seq_name]
    return(indexes)

def prepare_hto_index(index, path_quant_out):
    features_tsv = os.path.join(path_quant_out, 'features.tsv')
    t2g = os.path.join(path_quant_out, 't2g_hto.tsv')

    # The transcript2gene file for the HTO is a 1:1 mapping between the feature names
    with open(features_tsv, 'r') as in_file, open(t2g, 'w') as out_file:
        for line in in_file:
            feature_id = line.split('\t')[0]
            print(feature_id + '\t' + feature_id, file = out_file)

    logging.info('##\tbuilding HTO index')
    index_path = os.path.join(path_quant_out, 'hto_index')
    os.system(f'salmon index \
        -t {features_tsv} \
        -i {index_path} \
        --features \
        -k 7 \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')
    return index_path, t2g

def salmon_alevin(files_read_1, files_read_2, index_path, salmon_alevin_param, path_quant_out):
    read_1_string = ' '.join(files_read_1)
    read_2_string = ' '.join(files_read_2)
    logging.info(f'##\tmapping with alevin')
    os.system(f'salmon alevin \
        -l ISR \
        -i {index_path} \
        -1 {read_1_string} \
        -2 {read_2_string} \
        -p {snakemake.threads} \
        -o {path_quant_out}/map \
        --sketch \
        {salmon_alevin_param} \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')

def alevin_fry(path_quant_out, direction, tgmap):
    logging.info(f'##\tgenerating permit list')
    whitelist = snakemake.config['whitelist']
    os.system(f'alevin-fry generate-permit-list \
        -d {direction} \
        -i {path_quant_out}/map \
        -o {path_quant_out}/quant \
        -u {whitelist} \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')
    
    logging.info(f'##\tcollating')
    os.system(f'alevin-fry collate \
        -r {path_quant_out}/map \
        -i {path_quant_out}/quant \
        -t 1 \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')
    logging.info(f'##\tquantifying')
    os.system(f'alevin-fry quant \
        --threads {snakemake.threads} \
        --input-dir {path_quant_out}/quant \
        --output-dir {path_quant_out}/res \
        --resolution cr-like \
        --use-mtx \
        --tg-map {tgmap} \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')

def count_lines(file):
    i = -1
    with open(file) as f:
            for i, _ in enumerate(f):
                pass
    return i + 1

def process_index(index, seq_names, index_type):
    logging.info(f'##\t-\tproccessing index:\t{index}')
    path_quant_out = os.path.join(
        snakemake.config['out_path'], 
        snakemake.config['com_id'], 
        index_type, 
        index)

    seq_names_join = ', '.join(seq_names)
    logging.info(f'##\t\t(from {seq_names_join})')

    run_complete = os.path.isdir(os.path.join(path_quant_out, "res"))

    if run_complete:
        logging.info('##\t-\tpipestance exists - skipping')
        return

    os.system(f'mkdir -p {path_quant_out}')

    files_read_1 = list()
    files_read_2 = list()
    for seq_name in seq_names:
        flowID = seq_name.split('_')[3][1:]
        fastq_folder = os.path.join(
            snakemake.config['fastq_path'],
            seq_name,
            'fastq',
            flowID
            )
        files_read_1 += glob.glob(f'{fastq_folder}/*{index}*R1*.fastq.gz')
        files_read_2 += glob.glob(f'{fastq_folder}/*{index}*R2*.fastq.gz')
                
    if index_type == "hto":
        index_path, alevin_fry_tgmap = prepare_hto_index(index, path_quant_out)
        direction = "fw"
        salmon_alevin_param = f'\
            --read-geometry 2[1-15] \
            --bc-geometry 1[1-16] \
            --umi-geometry 1[17-26]'
        usa_mode = 'false'
    else:
        index_path = snakemake.config['salmon_index']
        direction = "both"
        alevin_fry_tgmap = snakemake.config['T3G']
        t2g = snakemake.config['T3G']
        salmon_alevin_param = f'\
            --chromiumV3 \
            --tgMap {t2g}'
        usa_mode = 'true'

    salmon_alevin(files_read_1, files_read_2, index_path, salmon_alevin_param, path_quant_out)
    alevin_fry(path_quant_out, direction, alevin_fry_tgmap)
    
    alevin_res = os.path.join(path_quant_out, 'res', 'alevin')
    count_genes = count_lines(alevin_res + '/quants_mat_cols.txt')
    count_cells = count_lines(alevin_res + '/quants_mat_rows.txt')

    with open(f'{path_quant_out}/res/meta_info.json', 'w') as f:
            f.write('{')
            f.write('  "alt_resolved_cell_numbers": [],')
            f.write('  "dump_eq": false,')
            f.write(f'  "num_genes": {count_genes},')
            f.write(f'  "num_quantified_cells": {count_cells},')
            f.write('  "resolution_strategy": "CellRangerLike",')
            f.write(f'  "usa_mode": {usa_mode}')
            f.write('}')


main()

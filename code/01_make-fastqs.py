# Imports
import os
from math import ceil
import logging
from bs4 import BeautifulSoup

def main():
    # Init
    logging.basicConfig(level=logging.INFO,
        filename=snakemake.log[0],
        format='%(message)s')
                    
    logging.getLogger().addHandler(logging.StreamHandler())

    print('#'*80)
    print('#####                                                                      #####')
    print('#####                      Demultiplexing FASTQ files                      #####')
    print('#####                                                       ** / *****     #####')
    print('#'*80)
    print('#####\n##')
    print('##\tDemultiplexing libraries:')

    bcl_folders = snakemake.config["pin_10x"].keys()
    if snakemake.config["pin_10x"] == "10x + HTO":
        hto_bcls = snakemake.config["pin_hto"].keys()
        bcl_folders = bcl_folders.union(hto_bcls)

    for bcl_folder in bcl_folders:
        print(f'##\t-\t{bcl_folder}')

        index_metrics = os.path.join(snakemake.config['fastq_path'], bcl_folder, 'interop', 'IndexMetricsOut.bin')
        if os.path.isfile(index_metrics):
            print('##\t-\tpipestance exists - skipping')
        else:
            bcl2fastq(bcl_folder)
        print('##\n###')
        print('#'*80)

    os.system(f'touch {snakemake.output[0]}')


def read_cycles(file):
    with open(file, 'r') as f:
        xml = f.read()
        soup = BeautifulSoup(xml, "xml")

        r1 = soup.find('Read1NumberOfCycles').text
        r2 = soup.find('Read2NumberOfCycles').text
        i1 = soup.find('IndexRead1NumberOfCycles').text
        i2 = soup.find('IndexRead2NumberOfCycles').text

    return {'r1' : r1, 'r2' : r2, 'i1' : i1, 'i2' : i2 }


def bcl2fastq(bcl_folder):
    # read/write needs less threads than processing
    rw_threads = ceil(snakemake.threads * 0.1)
    p_threads = ceil(snakemake.threads * 0.8)

    bcl_path = os.path.join(snakemake.config['bcl_path'], bcl_folder)
    fastq_path = os.path.join(snakemake.config['fastq_path'], bcl_folder, 'fastq')
    interop_path = os.path.join(snakemake.config['fastq_path'], bcl_folder, 'interop')
    sample_sheet = os.path.join(snakemake.config['fastq_path'], 'SampleSheet_' + bcl_folder + '.csv')
    log_path = os.path.join(snakemake.config['log_path'], snakemake.config['com_id'], bcl_folder + '_bcl2fastq.log')

    param_file = os.path.join(snakemake.config['bcl_path'], bcl_folder, 'RunParameters.xml')
    cycles = read_cycles(param_file)

    # In case we need to support other chemistries, make an if else statement to reflect the different parameters
    bases_mask = (
        'Y' + cycles['r1'] + ',' +
        'I' + cycles['i1'] + ',' +
        'I' + cycles['i2'] + ',' +
        'Y' + cycles['r2']
        )
    min_trim = 8
    min_adapt = 8
    
    os.system(f'bcl2fastq --use-bases-mask={bases_mask} \
                --create-fastq-for-index-reads \
                --minimum-trimmed-read-length={min_trim} \
                --mask-short-adapter-reads={min_adapt} \
                --ignore-missing-positions \
                --ignore-missing-controls \
                --ignore-missing-filter \
                --ignore-missing-bcls \
                -r {rw_threads} -p {p_threads} -w {rw_threads} \
                -R {bcl_path} \
                --output-dir={fastq_path} \
                --interop-dir={interop_path} \
                --sample-sheet={sample_sheet} \
                --barcode-mismatches=1 \
                --no-lane-splitting \
                2>&1 | tee {log_path} -a')

main()

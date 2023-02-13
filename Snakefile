rule all:
    input:
        'success/aggregate_complete.txt'

rule fastqs:
    input:
        'success/ready.txt'
    output:
        temp('success/make-fastqs_complete.txt')
    log:
        config['log_path'] + '/' + config['com_id'] + '/01_make-fastqs.log'
    threads:
        40
    script:
        'code/01_make-fastqs.py'

rule quant:
    input:
        rules.fastqs.output
    output:
        temp('success/map-quant_complete.txt')
    log:
        config['log_path'] + '/' + config['com_id'] + '/02_map-quant.log'
    threads:
        40
    script:
        'code/02_map-quant.py'

rule update:
    input:
        rules.quant.output
    output:
        temp('success/update-pool-table_complete.txt')
    log:
        config['log_path'] + '/' + config['com_id'] + '/03_update-pool-table.log'
    script:
        'code/03_update-pool-table.R'

rule seurat:
    input:
        rules.update.output
    output:
        temp('success/make-seurat_complete.txt')
    log:
        config['log_path'] + '/' + config['com_id'] + '/04_make-seurat.log'
    script:
        'code/04_make-seurat.R'

rule aggregate:
    input:
        rules.seurat.output
    output:
        temp('success/aggregate_complete.txt')
    log:
        config['log_path'] + '/' + config['com_id'] + '/05_aggregate.log'
    script:
        'code/05_aggregate.R'


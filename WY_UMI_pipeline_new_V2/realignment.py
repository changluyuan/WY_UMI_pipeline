from __future__ import print_function
import sys
import generate_slurm
import utils


# utils.makedir('bam_raw')

def realign_single(cfg, species='human'):
    genome = None
    if species == 'human':
        genome = '/gpfs/genomedb/b37/v37_decoy_plus_phage.fasta'
    elif species == 'mouse':
        genome = '/gpfs/genomedb/other_species/mouse/GRCh38.75_plus_phage.fasta'
    elif species == 'PCR':
        genome = '/home/changluyuan/gpfs/data/meth_MsssI_5mC/reference/PCR_reference.fa'
    with open(cfg) as cfg:
        first = cfg.readline()
    first = first.strip().split("\t")
    samplename = first[0]
    rg_id = first[1]
    jobname = 'realignment_{}'.format(rg_id)
    cmd = '/gpfs/bin/anaconda2-xiaojuan/bin/python /gpfs/bin/bwa-meth-0.2.0/bwameth.py --reference {} ' \
                  '-t 20 bam_mkdup/{}.mkdup.R1.fastq.gz bam_mkdup/{}.mkdup.R2.fastq.gz | ' \
'/gpfs/bin/samtools-1.3.1/samtools view -F 4 -h -b -@ 20 | /gpfs/bin/samtools-1.3.1/samtools sort -@ 20 -o bam_mkdup/{}.rmdup.bam && samtools index bam_mkdup/{}.rmdup.bam'.format(genome, rg_id, rg_id, rg_id, rg_id)
    return jobname, cmd

def realign_merge(cfg, species='human'):
    genome = None
    if species == 'human':
        genome = '/gpfs/genomedb/b37/v37_decoy_plus_phage.fasta'
    elif species == 'mouse':
        genome = '/gpfs/genomedb/other_species/mouse/GRCh38.75_plus_phage.fasta'
    elif species == 'PCR':
        genome = '/home/changluyuan/gpfs/data/meth_MsssI_5mC/reference/PCR_reference.fa'
    with open(cfg) as cfg:
        first = cfg.readline()
    first = first.strip().split("\t")
    samplename = first[0]
    rg_id = first[1]
    jobname = 'realignment_{}'.format(samplename)
    cmd = '/gpfs/bin/anaconda2-xiaojuan/bin/python /gpfs/bin/bwa-meth-0.2.0/bwameth.py --reference {} ' \
                  '-t 20 bam_mkdup/{}.merged.mkdup.R1.fastq.gz bam_mkdup/{}.merged.mkdup.R2.fastq.gz | ' \
'/gpfs/bin/samtools-1.3.1/samtools view -F 4 -h -b -@ 20 | /gpfs/bin/samtools-1.3.1/samtools sort -@ 20 -o bam_mkdup/{}.merged.rmdup.bam && samtools index bam_mkdup/{}.merged.rmdup.bam'.format(genome, samplename,samplename,samplename,samplename)
    return jobname, cmd


if __name__ == "__main__":
    print(align(sys.argv[1], species='human'))

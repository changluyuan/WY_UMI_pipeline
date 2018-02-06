#import glob
import utils
#import subprocess
import os

pwd = os.getcwd()

def extract_from_single(cfg, species):
    genome = None
    if species == 'human':
        genome = '/gpfs/genomedb/b37/v37_decoy_plus_phage.fasta'
    elif species == 'mouse':
        genome = '/gpfs/genomedb/other_species/mouse/GRCh38.75_plus_phage.fasta'
    elif species == 'PCR':
        genome = '/gpfs/users/changluyuan/data/meth_MsssI_5mC/reference/PCR_reference.fa'
    utils.makedir('CpG.methylKit')
    with open(cfg) as cfg:
        first = cfg.readline()
    first = first.strip().split("\t")
    rg_id = first[1]
    jobname = 'extract_{}'.format(rg_id)
    cmd = '/gpfs/bin/PileOMeth/PileOMeth extract -D 50000 --methylKit --CHH --CHG ' \
          '-o CpG.methylKit/{}.mkdup {} bam_mkdup/{}.rmdup.bam'.format(rg_id, genome, rg_id)
    return jobname, cmd


def extract_from_merge(cfg, species='human'):
    genome = None
    if species == 'human':
        genome = '/gpfs/genomedb/b37/v37_decoy_plus_phage.fasta'
    elif species == 'mouse':
        genome = '/gpfs/genomedb/other_species/mouse/GRCh38.75_plus_phage.fasta'
    elif species == 'PCR':
        genome = '/gpfs/users/changluyuan/data/meth_MsssI_5mC/reference/PCR_reference.fa'
    utils.makedir('CpG.methylKit')
    with open(cfg) as cfg:
        first = cfg.readline()
    first = first.strip().split("\t")
    samplename = first[0]
    jobname = 'extract_{}'.format(samplename)
    cmd = '/gpfs/bin/PileOMeth/PileOMeth extract -D 50000 --methylKit --CHH --CHG ' \
          '-o CpG.methylKit/{}.merged.mkdup {} bam_mkdup/{}.merged.rmdup.bam '.format(samplename, genome, samplename)
    return jobname, cmd


# def collect(cfg):
# for methfile in glob.glob('CpG.methylKit/*_CpG.methylKit'):
# subprocess.Popen(['ln', '-s', '{}'.format(methfile), '{}'.format(path)])

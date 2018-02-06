from __future__ import print_function
import sys
import generate_slurm
import utils


# utils.makedir('bam_raw')


def align(cfg, species='human'):
    genome = None
    if species == 'human':
        genome = '/gpfs/genomedb/b37/v37_decoy_plus_phage.fasta'
    elif species == 'mouse':
        genome = '/gpfs/genomedb/other_species/mouse/GRCh38.75_plus_phage.fasta'
    elif species == 'PCR':
        genome = '/home/changluyuan/gpfs/data/meth_MsssI_5mC/reference/PCR_reference.fa'
    cmd_list = []
    jobname_list = []
    with open(cfg) as cfg:
        for line in cfg.readlines():
            line = line.strip().split("\t")
            sample = line[0]
            rg_id = line[1]
            fq1 = line[2]
            fq2 = line[3]
            rg = '"@RG\\tID:{}\\tSM:{}\\tLB:{}"'.format(rg_id, sample, sample)
            cmd = '/gpfs/bin/anaconda2-xiaojuan/bin/python /gpfs/bin/bwa-meth-0.2.0/bwameth.py --reference {} ' \
                  '-t 20 --read-group {} fastq_trimmed/$(' \
                  'basename {} .fastq.gz).trimmed.fastq.gz fastq_trimmed/$(' \
                  'basename {} .fastq.gz).trimmed.fastq.gz | ' \
'/gpfs/bin/samtools-1.3.1/samtools view -h -b -@ 20 | /gpfs/bin/samtools-1.3.1/samtools sort -@ 20 -o {}.bam'.format(genome, rg, fq1, fq2, rg_id)
            jobname = 'map_{}'.format(rg_id)
            cmd_list.append(cmd)
            jobname_list.append(jobname)
        return jobname_list, cmd_list

            #s = generate_slurm.Slurm(jobname, {"ntasks-per-node": 20})
            #s.run(cmd, tries=1)

if __name__ == "__main__":
    print(align(sys.argv[1], species='human'))

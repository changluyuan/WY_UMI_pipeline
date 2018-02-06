import utils


def merge(cfg):
    utils.makedir('bam_merged')
    samples = list()
    bams = list()
    with open(cfg) as cfg:
        for line in cfg:
            line = line.strip().split("\t")
            samples.append(line[0])
            bams.append(line[1])
    samplename = samples[0]
    bamstring = '.bam '.join(bams)
    jobname = 'merge_{}'.format(samplename)
    cmd = 'sambamba merge -t 20 bam_merged/{}.merged.bam {}.bam'.format(
        samplename, bamstring)
    return jobname, cmd


def mkdup_from_single(cfg,minmem):
    utils.makedir('bam_raw')
    utils.makedir('bam_mkdup')
    with open(cfg) as cfg:
        first = cfg.readline()
    first = first.strip().split("\t")
    samplename = first[0]
    rg_id = first[1]
    jobname = 'mkdup_{}'.format(rg_id)
    cmd = 'mv {}.bam bam_raw && samtools index bam_raw/{}.bam && python /home/changluyuan/gpfs/bin/UMI_position_unique_bam_v7_new.py bam_raw/{}.bam bam_mkdup/{}.mkdup.bam {} && samtools sort bam_mkdup/{}.mkdup.bam > bam_mkdup/{}.mkdup.sort.bam && samtools index bam_mkdup/{}.mkdup.sort.bam && rm bam_mkdup/{}.mkdup.bam'.\
        format(rg_id, rg_id, rg_id, rg_id, minmem, rg_id,rg_id, rg_id, rg_id)
    return jobname, cmd


def mkdup_from_merge(cfg,minmem):
    utils.makedir('bam_mkdup')
    with open(cfg) as cfg:
        first = cfg.readline()
    first = first.strip().split("\t")
    samplename = first[0]
    jobname = 'mkdup_{}'.format(samplename)
    cmd = 'python /home/changluyuan/gpfs/bin/UMI_position_unique_bam_v7_new.py bam_merged/{}.merged.bam ' \
          'bam_mkdup/{}.merged.mkdup.bam {} && samtools sort bam_mkdup/{}.merged.mkdup.bam > bam_mkdup/{}.merged.mkdup.sort.bam && samtools index bam_mkdup/{}.merged.mkdup.sort.bam && rm bam_mkdup/{}.merged.mkdup.bam'.format(samplename, samplename, minmem, samplename,samplename, samplename, samplename)
    return jobname, cmd

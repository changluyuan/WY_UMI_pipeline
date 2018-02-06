import utils


def trim(cfg):
    utils.makedir('fastq_trimmed')
    cmd_list = []
    jobname_list = []
    with open(cfg) as cfg:
        for line in cfg.readlines():
            line = line.strip().split("\t")
            rg_id = line[1]
            fq1 = line[2]
            fq2 = line[3]
            jobname = 'trim_{}'.format(rg_id)
            cmd = 'cutadapt -e 0.1 -q 20 -O 1 -m 20 -a AGATCGGAAGAGC -A ' \
                  'GCGAATTTCGACG -U 14 -u 22 -o fastq_trimmed/$(basename {} ' \
                  '.fastq.gz).trimmed.fastq.gz -p fastq_trimmed/$(basename ' \
                  '{} .fastq.gz).trimmed.fastq.gz fastq_tag_header/{}.seq1.smi.fq.gz fastq_tag_header/{}.seq2.smi.fq.gz'.format(fq1,fq2,rg_id,rg_id)
            #s = generate_slurm.Slurm(jobname, {"ntasks-per-node": 1})
            #s.run(cmd, tries=1)
            cmd_list.append(cmd)
            jobname_list.append(jobname)
        return jobname_list, cmd_list




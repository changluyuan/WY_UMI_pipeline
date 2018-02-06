import utils


def tag_to_header(cfg):
    utils.makedir('fastq_tag_header')
    cmd_list = []
    jobname_list = []
    with open(cfg) as cfg:
        for line in cfg.readlines():
            line = line.strip().split("\t")
            rg_id = line[1]
            fq1 = line[2]
            fq2 = line[3]
            jobname1 = 'tag_to_header_{}'.format(rg_id)
            cmd1 = 'python ~/gpfs/bin/Duplex-Sequencing_WY/Nat_Protocols_Version/tag_to_header.py --infile1 {} --infile2 {} --taglen 8 --spacerlen 0 --outprefix fastq_tag_header/{}'.format(fq1, fq2, rg_id)
            #s = generate_slurm.Slurm(jobname, {"ntasks-per-node": 1})
            #s.run(cmd, tries=1)
            cmd_list.append(cmd1)
            jobname_list.append(jobname1)
        return jobname_list, cmd_list




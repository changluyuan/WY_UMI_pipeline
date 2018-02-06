import argparse
import glob
from collections import defaultdict

import trimming
import alignment
import markdup
import extractmeth
import utils
import generate_slurm
import tag_to_header
import realignment
parser = argparse.ArgumentParser(description='')
parser.add_argument('--species', )
parser.add_argument('--minmem',default=3)
parser.add_argument('config_file', )
parser.add_argument('--trim')
parser.add_argument('--align')
parser.add_argument('--merge')
parser.add_argument('--mkdup')
parser.add_argument('--extract')
parser.add_argument('--collect')
parser.add_argument('--tag_to_header')

args = parser.parse_args()

config_file = args.config_file
species = args.species
minmem = args.minmem

# split cfg file
utils.makedir('tmp')
sample_dict = defaultdict(list)
with open(config_file) as cfg:
    cfg = cfg.readlines()
    for line in cfg:
        line1 = line.strip().split("\t")
        sample_dict[line1[0]].append(line)

for key in sample_dict:
    with open('tmp/tmp.{}.cfg'.format(key), 'w') as output:
        for item in sample_dict[key]:
            output.write('{}'.format(item))


for i in glob.glob('tmp/tmp.*.cfg'):
    if utils.count_lines(i) > 1:
        step0jobid = []
        for row in range(0, utils.count_lines(i)):
            step0 = generate_slurm.Slurm(tag_to_header.tag_to_header(i)[0][row],
                                         {"ntasks-per-node": 1})
            step0jobid.append(step0.run(tag_to_header.tag_to_header(i)[1][row]))
        step1jobid = []
        for row in range(0, utils.count_lines(i)):
            step1 = generate_slurm.Slurm(trimming.trim(i)[0][row],
                                         {"ntasks-per-node": 1})
            step1jobid.append(step1.run(trimming.trim(i)[1][row],depends_on=step0jobid))
        step2jobid = []
        for row in range(0, utils.count_lines(i)):
            step2 = generate_slurm.Slurm(alignment.align(i, species)[0][row],
                                     {"ntasks-per-node": 20})
            step2jobid.append(step2.run(alignment.align(i, species)[1][row],
                                   depends_on=step1jobid))
        step3 = generate_slurm.Slurm(markdup.merge(i)[0], {"ntasks-per-node": 20})
        step3jobid = step3.run(markdup.merge(i)[1], depends_on=step2jobid)
        step4 = generate_slurm.Slurm(markdup.mkdup_from_merge(i,minmem)[0], {"ntasks-per-node": 20})
        step4jobid = step4.run(markdup.mkdup_from_merge(i,minmem)[1], depends_on= [step3jobid])
        step4_5 = generate_slurm.Slurm(realignment.realign_merge(i, species)[0], {"ntasks-per-node": 20})
        step4_5jobid = step4_5.run(realignment.realign_merge(i, species)[1],depends_on= [step4jobid])
        step5 = generate_slurm.Slurm(extractmeth.extract_from_merge(i, species)[
                                         0], {"ntasks-per-node": 1})
        step5jobid = step5.run(extractmeth.extract_from_merge(i, species)[1],
                               depends_on=[step4_5jobid])
    elif utils.count_lines(i) == 1:
        step0jobid = []
        step0 = generate_slurm.Slurm(tag_to_header.tag_to_header(i)[0][0],
                                         {"ntasks-per-node": 1})
        step0jobid.append(step0.run(tag_to_header.tag_to_header(i)[1][0]))
	step1jobid = []
        step1 = generate_slurm.Slurm(trimming.trim(i)[0][0],
                                     {"ntasks-per-node": 1}) #only one row
        step1jobid.append(step1.run(trimming.trim(i)[1][0],depends_on=step0jobid))
	step2jobid = []
        step2 = generate_slurm.Slurm(alignment.align(i, species)[0][0],
                                     {"ntasks-per-node": 20})
        step2jobid.append(step2.run(alignment.align(i, species)[1][0], depends_on=
            step1jobid))
        step3 = generate_slurm.Slurm(markdup.mkdup_from_single(i,minmem)[0],
                                     {"ntasks-per-node": 20})
        step3jobid = step3.run(markdup.mkdup_from_single(i,minmem)[1], depends_on=
            step2jobid)
        step3_4 = generate_slurm.Slurm(realignment.realign_single(i, species)[0], {"ntasks-per-node": 20})
        step3_4jobid = step3_4.run(realignment.realign_single(i, species)[1], depends_on=[step3jobid])
        step4 = generate_slurm.Slurm(extractmeth.extract_from_single(i,
                                                                     species)[0], {"ntasks-per-node": 1})
        step4jobid = step4.run(extractmeth.extract_from_single(i, species)[
                                   1], depends_on=[step3_4jobid])

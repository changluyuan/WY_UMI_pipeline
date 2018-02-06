from __future__ import print_function
import subprocess
import sys
import os
import utils


utils.makedir("log")


TMPL_long = """\
#!/bin/bash

#SBATCH -J {name}
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -o log/{name}.%j.out
#SBATCH -e log/{name}.%j.err
#SBATCH --no-requeue
#SBATCH -A cnl
{header}
set -eo pipefail -o nounset
cd {workdir}

__script__"""


TMPL_short = """\
#!/bin/bash

#SBATCH -J {name}
#SBATCH -p cn-medium
#SBATCH -N 1
#SBATCH -o log/{name}.%j.out
#SBATCH -e log/{name}.%j.err
#SBATCH --no-requeue
#SBATCH -A cnm
{header}
set -eo pipefail -o nounset
cd {workdir}

__script__"""


class Slurm(object):
    def __init__(self, name, slurm_kwargs=None, tmpl=None,
                 scripts_dir="jobs", workdir=os.getcwd()):
        if slurm_kwargs is None:
            slurm_kwargs = {}
        if tmpl is None:
            tmpl = TMPL_short
        else:
            tmpl = TMPL_long

        header = []
        for k, v in slurm_kwargs.items():
            if len(k) > 1:
                k = "--" + k + "="
            else:
                k = "-" + k + " "
            header.append("#SBATCH %s%s" % (k, v))
        self.header = "\n".join(header)
        self.name = name
        self.tmpl = tmpl
        self.slurm_kwargs = slurm_kwargs
        utils.makedir("jobs")
        self.scripts_dir = os.path.abspath(scripts_dir)
        self.workdir = os.getcwd()

    def __str__(self):
        return self.tmpl.format(name=self.name, header=self.header,
                                workdir=self.workdir)

    def _tmpfile(self):
        return '{}/{}.job'.format(self.scripts_dir, self.name)

    def run(self, command, cmd_kwargs=None, _cmd="sbatch", depends_on=None):

        if cmd_kwargs is None:
            cmd_kwargs = {}

        args = []
        for k, v in cmd_kwargs.items():
            args.append("export %s=%s" % (k, v))
        args = "\n".join(args)

        tmpl = str(self).replace("__script__", args + "\n###\n" + command +
                                 "\n")
        if depends_on is None or (len(depends_on) == 1 and depends_on[0] is None):
            depends_on = []

        with open(self._tmpfile(), "w") as sh:
            sh.write(tmpl)

        # job_id = None
        args = [_cmd]
        args.extend([("--dependency=afterok:%d" % int(d)) for d in depends_on])
        args.append(sh.name)
        res = subprocess.check_output(args).strip()
        print(res, file=sys.stderr)
        if not res.startswith(b"Submitted batch"):
            return None
        job_id = int(res.split()[-1])
        return job_id

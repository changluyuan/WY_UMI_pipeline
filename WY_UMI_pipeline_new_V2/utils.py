import os
import subprocess


def makedir(path):
    path = path.strip()
    is_exits = os.path.exists(path)
    if not is_exits:
        os.makedirs(path)


def count_lines(filename):
    out = subprocess.Popen(['wc', '-l', filename], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT).communicate()[0]
    return int(out.partition(b' ')[0])

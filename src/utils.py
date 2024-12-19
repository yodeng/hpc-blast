#!/usr/bin/env python

import os
import re
import sys
import gzip
import shlex
import shutil
import signal
import tempfile
import argparse
import subprocess

from runjob.parser import *
from runjob import log, runsge
from runjob.config import Config

from ._version import __version__


class Zopen(object):

    def __init__(self, name,  mode="rb"):
        self.name = name
        self.mode = mode
        self.handler = None

    def __enter__(self):
        if self.name.endswith(".gz"):
            if "r" in self.mode:
                p = subprocess.Popen(
                    ["gzip", "-c", "-d", self.name], stdout=subprocess.PIPE)
                self.handler = p.stdout
            else:
                self.handler = gzip.open(self.name, self.mode)
        else:
            self.handler = open(self.name, self.mode)
        return self.handler

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.handler:
            self.handler.close()


class MultiFileOpen(object):

    def __init__(self, mode="rb", **infiles):
        self.info = infiles
        self.handler = {}
        self.mode = mode

    def __enter__(self):
        for k, f in self.info.items():
            if f.endswith(".gz"):
                self.handler[int(k)] = gzip.open(f, self.mode)
            else:
                self.handler[int(k)] = open(f, self.mode)
        return self.handler

    def __exit__(self, exc_type, exc_val, exc_tb):
        for _, h in self.handler.items():
            h.close()


class ArgumentsError(Exception):
    pass


def mkdir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def get_fastx_type(fastx):
    fastx_ = fastx.lower()
    if fastx_.endswith(".gz"):
        fastx_ = fastx_[:-3]
    if fastx_.endswith(".fasta") or fastx_.endswith(".fa"):
        return "fasta"
    elif fastx_.endswith(".fastq") or fastx_.endswith(".fq"):
        return "fastq"
    fastx_type = "fastq"
    with Zopen(fastx) as fi:
        curr_line = 0
        for line in fi:
            if line.startswith(b">"):
                return "fasta"
            elif curr_line == 4:
                return "fastq"
            curr_line += 1
        if curr_line == 0:
            sys.exit("No entries in %s" % (fastx))
    return fastx_type


def which(program, paths=None):
    ex = os.path.dirname(sys.executable)
    found_path = None
    fpath, fname = os.path.split(program)
    if fpath:
        program = canonicalize(program)
        if is_exe(program):
            found_path = program
    else:
        if is_exe(os.path.join(ex, program)):
            return os.path.join(ex, program)
        paths_to_search = []
        if isinstance(paths, (tuple, list)):
            paths_to_search.extend(paths)
        else:
            env_paths = os.environ.get("PATH", "").split(os.pathsep)
            paths_to_search.extend(env_paths)
        for path in paths_to_search:
            exe_file = os.path.join(canonicalize(path), program)
            if is_exe(exe_file):
                found_path = exe_file
                break
    return found_path


def is_exe(file_path):
    return (
        os.path.exists(file_path)
        and os.access(file_path, os.X_OK)
        and os.path.isfile(os.path.realpath(file_path))
    )


def human_size_parse(size):
    s, u = re.search("(\d+(?:\.\d+)?)(\D*)", str(size)).group(1, 2)
    s = float(s)
    if s < 1 and not u:
        u = "M"
    if u:
        for unit in ['B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z']:
            if u.upper().strip()[0] == unit:
                return int(s)
            s *= 1024
    else:
        return int(s)


def callcmd(cmd, run=True, verbose=False):
    if not cmd:
        return
    if not run:
        if verbose:
            print(cmd)
        return
    if verbose:
        print(cmd)
        subprocess.check_call(cmd, shell=True, stdout=sys.stdout,
                              stderr=sys.stderr)
    else:
        with open(os.devnull, "w") as fo:
            subprocess.check_call(cmd, shell=True, stdout=fo, stderr=fo)


def canonicalize(path):
    return os.path.abspath(os.path.expanduser(path))


blast_dbtype = {
    "blastn": "nucl",
    "blastp": "prot",
    "blastx": "prot",
    "tblastn": "nucl",
    "tblastx": "nucl",
}


def rate_parser(parser):
    rate_args = parser.add_argument_group("rate arguments")
    rate_args.add_argument('--retry', help="retry N times of the error job, 0 or minus means do not re-submit.",
                           type=int, default=0, metavar="<int>")
    rate_args.add_argument('--retry-sec', help="retry the error job after N seconds.",
                           type=int, default=2, metavar="<int>")
    rate_args.add_argument('--max-check', help="maximal number of job status checks per second, fractions allowed.",
                           type=float, default=DEFAULT_MAX_CHECK_PER_SEC, metavar="<float>")
    rate_args.add_argument('--max-submit', help="maximal number of jobs submited per second, fractions allowed.",
                           type=float, default=DEFAULT_MAX_SUBMIT_PER_SEC, metavar="<float>")


def resource_parser(parser):
    resource_args = parser.add_argument_group("resource arguments")
    resource_args.add_argument("--queue", type=str, help="queue/partition for running, multi-queue can be sepreated by whitespace. (default: all accessed)",
                               nargs="*", metavar="<queue>")
    resource_args.add_argument("--node", type=str, help="node for running, multi-node can be sepreated by whitespace. (default: all accessed)",
                               nargs="*", metavar="<node>")
    resource_args.add_argument("--round-node", action="store_true", help="round all define node per job for load balance",
                               default=False)
    resource_args.add_argument("--cpu", type=int,
                               help="max cpu number used.", default=1, metavar="<int>")
    resource_args.add_argument("--memory", type=int,
                               help="max memory used (GB).", default=1, metavar="<int>")


def control_parser(parser):
    control_args_parser = parser.add_argument_group("control arguments")
    ex_args = control_args_parser.add_mutually_exclusive_group(required=False)
    ex_args.add_argument("--split", type=int, default=10,
                         help='split query into num of chunks, 10 by default', metavar="<int>")
    ex_args.add_argument("--size", type=int,
                         help='split query into multi chunks with N sequences', metavar="<int>")
    ex_args.add_argument("--filesize", type=str,
                         help="split query into multi chunks with define filesize, 1G, 500M", metavar="<str/float>")
    control_args_parser.add_argument("--num", type=int,
                                     help='max number of chunks run parallelly, all chunks by default', metavar="<int>")


def hpcblast_rate_resource_args(args):
    parser = argparse.ArgumentParser(allow_abbrev=False)
    resource_parser(parser)
    rate_parser(parser)
    # control_parser(parser)
    out = []
    args_map = args.__dict__ if isinstance(
        args, argparse.Namespace) else dict(args)
    for act in parser._actions:
        dest = act.dest
        if dest in args_map and args_map[dest] != act.default:
            value = args_map[dest]
            if isinstance(value, list):
                value = " ".join(value)
            dest = dest.replace("_", "-")
            if act.nargs == 0:
                out.append(f"--{dest}")
            else:
                out.extend([f"--{dest}", value])
    return " ".join(map(str, out))


def HPCBlastArg():
    parser = argparse.ArgumentParser(
        description="hpc-blast <OPTIONS> <blast command>",
        formatter_class=CustomHelpFormatter,
        add_help=False, allow_abbrev=False)
    control_parser(parser)
    parser.add_argument("--tempdir", type=str, required=False,
                        help='hpc blast temp directory', metavar="<dir>")
    parser.add_argument("--log", type=str,
                        help='append hpc-blast log info to file, sys.stdout by default', metavar="<file>")
    parser.add_argument("--local", action='store_true',
                        help="run blast in localhost instead of sge", default=False)
    parser.add_argument("--slurm", action='store_true',
                        help="run blast in slurm", default=False)
    parser.add_argument('--debug', action='store_true',
                        help='log debug', default=False)
    parser.add_argument('--version',
                        action='version', version="v" + __version__)
    parser.add_argument("-h", '--help',
                        action='help', help="show this help message and exit")
    parser.add_argument("blast", type=str,
                        help='blast command, required', metavar="<blast command>")
    resource_parser(parser)
    rate_parser(parser)
    args, unknown_args = parser.parse_known_args()
    args.blast_db = []
    c, o, d, q = 0, 0, 0, 0
    db_end = False
    for a in sys.argv[1:]:
        if d and a.startswith("-"):
            db_end = True
            d = 0
        if a == "-num_threads":
            c = 1
        elif a == "-out":
            unknown_args.remove(a)
            o = 1
        elif a == "-db":
            unknown_args.remove(a)
            d = 1
        elif a == "-query":
            unknown_args.remove(a)
            q = 1
        elif c:
            args.cpu = args.cpu or int(a)
            c = 0
        elif o:
            unknown_args.remove(a)
            args.outfile = a
            o = 0
        elif d:
            args.blast_db.append(a)
            unknown_args.remove(a)
        elif q:
            unknown_args.remove(a)
            args.query = a
            q = 0
        elif hasattr(args, "cpu") and hasattr(args, "outfile") and hasattr(args, "blast_db") and hasattr(args, "query") and db_end:
            break
    if not os.path.basename(args.blast).startswith("blast"):
        parser.error("hpc-blast argument error")
    return args, unknown_args

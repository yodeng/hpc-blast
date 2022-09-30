#!/usr/bin/env python

import os
import sys
import gzip
import shlex
import argparse
import subprocess

from random import randrange

from runjob.utils import Mylog
from runjob.config import Config
from runjob.sge_run import RunSge
from runjob.sge import ParseSingal

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


def mkdir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


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


def callcmd(cmd, run=True, verbose=False):
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


blast_dbtype = {"blastn": "nucl",
                "blastp": "prot",
                "blastx": "prot",
                "tblastn": "nucl",
                "tblastx": "nucl",
                }


def HPCBlastArg():
    parser = argparse.ArgumentParser(
        description="hpc-blast <OPTIONS> <blast options>")
    parser.add_argument("--blast", type=str, required=True,
                        help='blast type or blast bin path', metavar="<str>")
    parser.add_argument("--split", type=int, default=10,
                        help='split query into num of chunks, 10 by default', metavar="<int>")
    parser.add_argument("--queue", type=str, default="all.q",
                        help='sge queue, all.q by default, multi-queue can be sepreated by whitespace', nargs="*", metavar="<str>")
    parser.add_argument("--cpu", type=int, default=1,
                        help='cpu usage for sge, 1 by default', metavar="<int>")
    parser.add_argument("--memory", type=int, default=1,
                        help='memory (GB) usage for sge, 1 by default', metavar="<int>")
    parser.add_argument("--num", type=int,
                        help='max number of chunks run parallelly, default: all', metavar="<int>")
    parser.add_argument("--output", type=str, required=True,
                        help='hpc blast output directory', metavar="<str>")
    parser.add_argument("--log", type=str,
                        help='append hpc-blast log info to file, sys.stdout by default', metavar="<file>")
    parser.add_argument("--local", action='store_true',
                        help="run blast in localhost instead of sge", default=False)
    parser.add_argument('--version',
                        action='version', version="v" + __version__)
    args, unknown_args = parser.parse_known_args()
    c, o, d, q = 0, 0, 0, 0
    for a in sys.argv[1:]:
        if a == "-num_threads":
            c = 1
        elif a == "-out":
            unknown_args.remove(a)
            o = 1
        elif a == "-db":
            d = 1
        elif a == "-query":
            unknown_args.remove(a)
            q = 1
        elif c:
            if not args.cpu:
                args.cpu = int(a)
            c = 0
        elif o:
            unknown_args.remove(a)
            args.outfile = a
            o = 0
        elif d:
            args.blast_db = a
            d = 0
        elif q:
            unknown_args.remove(a)
            args.query = a
            q = 0
        elif hasattr(args, "cpu") and hasattr(args, "outfile") and hasattr(args, "blast_db") and hasattr(args, "query"):
            break
    return args, unknown_args

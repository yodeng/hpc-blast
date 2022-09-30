#!/usr/bin/env python

from .utils import *


class HPCBlast(object):
    def __init__(self, args=None, blast_options=None):
        self.args = args
        self.outdir = args.output
        self.args.mode = "sge"
        if self.args.local:
            self.args.mode = "local"
        self.args.jobname = "hpc_blast_%d" % os.getpid()
        self.args.logdir = os.path.join(self.outdir, "logs")
        self.args.workdir = os.getcwd()
        self.args.startline = 0
        self.args.groups = 1
        self.btype = args.blast
        self.db = args.blast_db
        self.blast_exe = which(self.btype)
        self.query = args.query
        self.outfile = args.outfile
        self.blast_options = blast_options
        self.chunk_files = {}
        self.chunk_res = []
        self.blast_scripts = ""

    def split_fasta(self, num=10):
        self.chunk_files = {str(i): os.path.join(
            self.outdir, "chunks", "split.%05d.fa" % i) for i in range(num)}
        if os.path.isfile(self.query):
            mkdir(os.path.join(self.outdir, "chunks"))
        with Zopen(self.query, mode="rb") as fi, MultiFileOpen(mode="wb", **self.chunk_files) as fo:
            for line in fi:
                if line.startswith(b">"):
                    n = randrange(0, num)
                    fh = fo[n]
                fh.write(line)

    def cache_blast_db_cmd(self):
        blastdb_path = which("blastdb_path")
        vmtouch = which("vmtouch")
        parallel = which("parallel")
        cmd = "%s -db %s -dbtype %s -getvolumespath | tr ' ' '\n' | %s %s -tqm 5G 2> /dev/null" % (
            blastdb_path, self.db, blast_dbtype[os.path.basename(self.btype)], parallel, vmtouch)
        return cmd

    def write_blast_sh(self, out="hpc_blast.sh"):
        self.blast_scripts = os.path.join(self.outdir, out)
        with open(self.blast_scripts, "w") as fo:
            for _, fa in self.chunk_files.items():
                if not os.path.getsize(fa):
                    continue  # ignore empty query
                name = os.path.basename(fa).split(".")
                result = os.path.join(
                    self.outdir, "results", "result."+name[1])
                self.chunk_res.append(result)
                cmdline = [self.blast_exe, ] + self.blast_options
                cmdline.extend(["-out", result, "-query", fa])
                fo.write(shlex.join(cmdline)+"\n")

    def run_blast(self):
        mkdir(os.path.join(self.outdir, "results"))
        self.args.jobfile = self.blast_scripts
        self.args.force = True  # force to resubmit
        conf = Config()
        conf.update_dict(**self.args.__dict__)
        if os.path.isfile(self.blast_scripts):
            mkdir(os.path.join(self.outdir, "logs"))
            runsge_loger = Mylog(self.args.log, "info", name=RunSge.__module__)
            runsge = RunSge(config=conf)
            h = ParseSingal(obj=runsge, name=self.args.jobname,
                            mode=self.args.mode, conf=conf)
            h.start()
            runsge.run(times=0)
            if not runsge.sumstatus():
                os.kill(os.getpid(), signal.SIGTERM)

    def mergs_res(self):
        if self.chunk_res:
            mergs = "cat %s > %s" % (" ".join(self.chunk_res), self.outfile)
            callcmd(mergs)

    def run(self):
        splits = self.args.split
        self.split_fasta(splits)
        self.write_blast_sh()
        self.run_blast()
        self.mergs_res()

#!/usr/bin/env python

from .utils import *

__all__ = ["HPCBlast", "HPCBlastArg"]


class HPCBlast(object):

    def __init__(self, args=None, blast_options=None):
        self.outfile = os.path.abspath(args.outfile)
        self.tempdir = os.path.abspath(args.tempdir or tempfile.mktemp(
            prefix="hpc-blast_", dir=os.path.dirname(self.outfile)))
        self.args = args
        self._create_sge_args()
        self.btype = args.blast
        self.db = args.blast_db
        self.blast_exe = which(self.btype)
        self.query = args.query
        self.blast_options = blast_options
        self.chunk_files = {}
        self.chunk_res = []
        self.blast_scripts = ""
        self.finished = False
        if not self.blast_exe or not "blast" in os.path.basename(self.blast_exe):
            raise ArgumentsError(
                "blast not found in this environment")
        self.cleandir = not args.tempdir

    def _quotation_outfmt(self):
        if "-outfmt" in self.blast_options:
            fmt_idx = self.blast_options.index("-outfmt")
            op = []
            s = fmt_idx+1
            e = len(self.blast_options)
            if len(self.blast_options[s].split()) == 1:
                for n, i in enumerate(self.blast_options[s:]):
                    if i.startswith("-"):
                        e = s + n - 1
                        break
                    op.append(i.strip("'").strip('"'))
                else:
                    e = s + n
                self.blast_options[s:e+1] = [" ".join(op)]

    def _create_sge_args(self):
        self.args.mode = "sge"
        if self.args.local:
            self.args.mode = "local"
        self.args.logdir = os.path.join(self.tempdir, "logs")
        self.args.workdir = os.getcwd()
        self.args.startline = 0
        self.args.groups = 1

    def split_fastx_by_size(self, size=0):
        self.chunk_files = {}
        if size <= 0:
            return
        if os.path.isfile(self.query):
            mkdir(os.path.join(self.tempdir, "chunks"))
        fx = get_fastx_type(self.query)
        s = fh = 0
        with Zopen(self.query, mode="rb") as fi:
            if fx == "fasta":
                for line in fi:
                    if line.startswith(b">"):
                        n, mod = divmod(s, size)
                        s += 1
                        if mod == 0:
                            if fh:
                                fh.close()
                            fo = os.path.join(
                                self.tempdir, "chunks", "split.%05d.fa" % n)
                            self.chunk_files[n] = fo
                            fh = open(fo, "wb")
                    fh.write(line)
            elif fx == "fastq":
                for i, line in enumerate(fi):
                    x = i % 4
                    if x == 0:
                        n, mod = divmod(s, size)
                        s += 1
                        if mod == 0:
                            if fh:
                                fh.close()
                            fo = os.path.join(
                                self.tempdir, "chunks", "split.%05d.fa" % n)
                            self.chunk_files[n] = fo
                            fh = open(fo, "wb")
                        line = b">" + line[1:]
                    elif x > 1:
                        continue
                    fh.write(line)
        if fh:
            fh.close()

    def split_fastx_by_part(self, part=10):
        self.chunk_files = {str(i): os.path.join(
            self.tempdir, "chunks", "split.%05d.fa" % i) for i in range(part)}
        if os.path.isfile(self.query):
            mkdir(os.path.join(self.tempdir, "chunks"))
        fx = get_fastx_type(self.query)
        num = 0
        with Zopen(self.query, mode="rb") as fi, MultiFileOpen(mode="wb", **self.chunk_files) as fo:
            if fx == "fasta":
                for line in fi:
                    if line.startswith(b">"):
                        fh = fo[num % part]
                        num += 1
                    fh.write(line)
            elif fx == "fastq":
                for i, line in enumerate(fi):
                    x = i % 4
                    if x == 0:
                        fh = fo[num % part]
                        num += 1
                        line = b">" + line[1:]
                    elif x > 1:
                        continue
                    fh.write(line)

    @property
    def cache_blast_db(self):
        blastdb_path = which("blastdb_path")
        vmtouch = which("vmtouch")
        cmd = ""
        if vmtouch and blastdb_path:
            cmd = "%s -db %s -dbtype %s -getvolumespath | xargs %s -tqm 5G 2> /dev/null" % (
                blastdb_path, self.db, blast_dbtype[os.path.basename(self.btype)], vmtouch)
        return cmd

    def write_blast_sh(self, out="hpc_blast.sh"):
        self.blast_scripts = os.path.join(self.tempdir, out)
        self._quotation_outfmt()
        with open(self.blast_scripts, "w") as fo:
            for _, fa in self.chunk_files.items():
                if not os.path.getsize(fa):
                    continue  # ignore empty query file
                name = os.path.basename(fa).split(".")
                result = os.path.join(
                    self.tempdir, "results", "result."+name[1])
                self.chunk_res.append(result)
                cmdline = [self.blast_exe, ] + self.blast_options
                cmdline.extend(["-out", result, "-query", fa])
                fo.write(shlex.join(cmdline)+"\n")

    def run_blast(self):
        mkdir(os.path.join(self.tempdir, "results"))
        self.args.jobfile = self.blast_scripts
        self.args.force = True  # force to resubmit
        conf = Config()
        conf.update_dict(**self.args.__dict__)
        if os.path.isfile(self.blast_scripts):
            mkdir(os.path.join(self.tempdir, "logs"))
            loger = log(self.args.log, "info")
            job = runsge(config=conf)
            job.set_rate(20)
            job.run()
            loger.info("hpc blast finished")

    def mergs_res(self, block_size=65536):
        if self.chunk_res:
            with open(self.outfile, "wb") as fo:
                for f in self.chunk_res:
                    with open(f, "rb") as fi:
                        buf = fi.read(block_size)
                        while buf:
                            fo.write(buf)
                            buf = fi.read(block_size)
            self.finished = True

    def run(self):
        if self.args.size:
            self.split_fastx_by_size(self.args.size)
        else:
            self.split_fastx_by_part(self.args.split)
        self.write_blast_sh()
        self.run_blast()
        self.mergs_res()

    def __del__(self):
        tempath = ["chunks", "results", "logs", "hpc_blast.sh"]
        try:
            if self.finished and os.path.isdir(self.tempdir):
                for p in tempath:
                    path = os.path.join(self.tempdir, p)
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                else:
                    if not os.listdir(self.tempdir):
                        shutil.rmtree(self.tempdir)
        except:
            pass

# hpcblast

`hpcblast` software provides a simple and efficient method for running the `NCBI-BLAST+` suite program in localhost or the `HPC` environment (Sun Grid Engine). It splits the input sequence file and run all chunked tasks in parallel to accelerate the blast search speed. When there are many sequences in the input file for blast comparison and the running speed is slow, using hpcblast can significantly improve the performance. All of this can be easy done only by adding the `hpc-blast` command at the head of your blast command line.



### Features

+ hpcblast splits the input sequence file into small files and runs all tasks in parallel.
+ hpcblast supports fasta/fastq sequence format file input and gzip compression allowed, there is no need to decompress fastq and convert it to fasta for blast .
+ hpcblast manages and schedules all tasks by [**runjob**](https://github.com/yodeng/runjob).
+ hpcblast is compatible with all `NCBI-BLAST+` options, and all results are the same except the order of output aligned segment.
+ the `hpc-blast` option uses **two** `-` flag, while `NCBI-BLAST+` options use **one** `-` flag.
+ you need to pay attention to the CPU and memory resources used In the parallelized environment.
+ other efficient  features ...



### Requirements

+ python >= 3.5
+ [ncbi-blast+](https://anaconda.org/bioconda/blast)



### Installation

> from git repo (for recommend)

```
pip3 install git+https://github.com/yodeng/hpc-blast.git
```

> from pypi:

```
pip3 install hpcblast -U
```



### Usage

```
$ hpc-blast --help 
usage: hpc-blast [--split <int> | --size <int>] [--num <int>] [--tempdir <dir>] [--log <file>] [--local]
                 [--version] [-h] [--queue [<str> ...]] [--cpu <int>] [--memory <int>]
                 <blast command>

hpc-blast <OPTIONS> <blast command>

positional arguments:
  <blast command>      blast command, required

optional arguments:
  --split <int>        split query into num of chunks, 10 by default
  --size <int>         split query into multi chunks with N sequences
  --num <int>          max number of chunks run in parallel, all chunks by default
  --tempdir <dir>      hpc blast temp directory
  --log <file>         append hpc-blast log info to file, sys.stdout by default
  --local              run blast in localhost instead of sge
  --version            show program's version number and exit
  -h, --help           show this help message and exit

sge arguments:
  --queue [<str> ...]  sge queue, multi-queue can be sepreated by whitespace, all access queue by default
  --cpu <int>          cpu usage for sge, 1 by default, max(--cpu, -num_threads) will be used
  --memory <int>       memory (GB) usage for sge, 1 by default
```



### Example

> in Sun Grid Engine:

```
hpc-blast blastn -query test.fastq.gz -db /data/refdb -num_threads 4 -outfmt 6 qseqid qlen qstart qend -out test.m6
```

> in localhost:

```
hpc-blast --local blastn -query test.fastq.gz -db /data/refdb -num_threads 4 -outfmt 6 qseqid qlen qstart qend -out test.m6
```


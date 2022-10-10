# hpcblast

`hpcblast`软件用于在`HPC`环境下运行`NCBI-BLAST+`套件程序，提供了简单高效的方法，通过将比对输入文件进行拆分并利用`HPC`环境加速blast查找速度。当`blast`比对输入序列较多，运行速度较慢时使用`hpcblast`有明显的性能提升。

### 1. 特性

+ `hpcblast`将比对序列随机拆分成小文件并将任务并发运行在`SGE`或本地环境中
+ `hpcblast`兼容全部的`NCBI-BLAST+`选项和命令，结果和`NCBI-BLAST+`相同
+ `hpcblast`支持`fa/fq`格式，支持`gz`压缩格式输入，无需将`fastq`解压并转为`fasta`步骤
+ `hpcblast`统一的任务管理调度，自动管理资源，当程序中断或报错情况，自动销毁运行中的任务
+ `hpcblast`支持自定义并发控制和资源申请

### 2. 依赖

##### 2.1 运行环境

+ linux64
+ python >= 3.8
+ [ncbi-blast+](https://anaconda.org/bioconda/blast)

##### 2.2 其他python模块依赖

+ [Cython](https://github.com/cython/cython)
+ [runjob](https://github.com/yodeng/runjob)

### 3. 安装

> git仓库安装 (for recommend)

```
pip3 install git+https://github.com/yodeng/hpc-blast.git
```

> Pypi安装

```
pip3 install hpcblast -U
```

### 4. 使用

```
$ hpc-blast -h 
usage: hpc-blast [--split <int>] [--queue [<str> ...]] [--cpu <int>] [--memory <int>] [--num <int>] --output <str> [--log <file>] [--local] [--version] [-h] <blast command>

hpc-blast <OPTIONS> <blast commands>

positional arguments:
  <blast command>      blast command

optional arguments:
  --split <int>        split query into num of chunks, 10 by default
  --queue [<str> ...]  sge queue, all.q by default, multi-queue can be sepreated by whitespace
  --cpu <int>          cpu usage for sge, 1 by default, max(--cpu, -num_threads) will be used
  --memory <int>       memory (GB) usage for sge, 1 by default
  --num <int>          max number of chunks run parallelly, default: all
  --output <str>       hpc blast output directory
  --log <file>         append hpc-blast log info to file, sys.stdout by default
  --local              run blast in localhost instead of sge
  --version            show program's version number and exit
  -h, --help           show this help message and exit
```

相关参数解释如下：

| 参数       | 描述                                                         |
| ---------- | ------------------------------------------------------------ |
| --split    | 将输入序列拆分为多少份                                       |
| --queue    | 如果运行在`sge`环境，则传入指定的队列名，多个队列用空白隔开，默认`all.q` |
| --cpu      | 如果运行在`sge`环境，则任务申请的`cpu`资源数，默认1，会识别`blast`参数中的`-num_thread`参数，使用较大的值 |
| --memory   | 如果运行在`sge`环境，则任务申请的内存，默认1，单位`GB`       |
| --num      | `hpcblast`任务队列，将拆分的块发送到`hpcblast`任务队列中，代表同时运行的拆分块任务数，默认全部拆分块同时运行 |
| --output   | `hpcblast`输出目录，不存在会自动创建                         |
| --log      | `hpcblast`输出日志文件，默认标准输出                         |
| --local    | 强制`hpcblast`本地并发运行，不使用`sge`运行。                |
| --version  | 打印`hpcblast`版本并退出                                     |
| -h, --help | 打印`hpcblast`帮助并退出                                     |

**说明**：

+ `hpc-blast`选项使用两个中划线`--`，`blast`选项使用一个中划线`-`
+ `hpc-blast`命令和选项后面直接跟需要运行的blast命令即可，完全兼容，命令阻塞直到全部任务运行完成
+ `hpcblast`任务队列为最大允许同时运行的任务块，当某个任务块运行完成时会立即添加一个任务块到任务队列中并启动运行
+ `hpcblast`主线程管理各并发子任务块，主线程几乎不占用cpu和内存资源，主程序终止，其管理的本地或sge任务自动清空，无需手动销毁
+ `hpc-blast`识别`blast`命令中的`-query`和`-num_thread`参数，`-query`后面可直接跟`fa/fq`格式文件（支持`*.gz`压缩格式）
+ `hpc-blast`随机拆分输入文件，拆分过程无须将全部序列读入内存，无内存消耗，比对输出结果无序
+ 并发环境，需注意使用的`cpu`和内存资源

### 5. 示例

以下命令表示将`test.fastq.gz`序列随机拆分成`200`份进行`blastn`比对，任务投递到`sge`集群中，投递队列为`all.q`和`test.q`，使用资源`virtual_free=1g,num_proc=4`，允许同时运行的最大任务数为`50`个，`blastn`输出结果文件为`test.blast.m6`，`hpc-blast`输出目录为当前目录下的`out`目录。

```
hpc-blast --split 200 --queue all.q test.q --num 50 --output out blastn -query test.fastq.gz -db /data/refdb -num_threads 4 -out test.blast.m6 -outfmt "6 qseqid qlen qstart qend"
```

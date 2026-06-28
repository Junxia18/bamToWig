/*
bamToWig $2014/12/09/$ @Jian-Hua Yang yangjh7@mail.sysu.edu.cn
*/
#include "BamAux.h"
#include "BamReader.h"
#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace BamTools;
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <locale>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#include "bioUtils.h"
#include "faiFile.h"
#include "bedFile.h"
#include "varFile.h"
#include "isoFile.h"
#include "samFile.h"
#include "homer_statistics.h"
#include "statistic.h"

#include "bamToWig.h"

/* [FIX-B] 引入 POSIX 目录操作头文件，替代存在安全隐患的 system("mkdir -p ...")
 */
#include <sys/stat.h>
#include <sys/types.h>

/* [新增] 跨平台递归建目录函数，相当于 "mkdir -p"。
 * 需要 Linux POSIX mkdir 接口（GCC 4.8.5 环境完全支持）。
 * 返回值：0 成功，-1 失败。*/
static int mkdirRecursive(const std::string& path) {
  if (path.empty())
    return 0;
  struct stat st;
  /* 如果目录已存在，直接返回成功 */
  if (stat(path.c_str(), &st) == 0)
    return 0;
  /* 先递归创建父目录 */
  size_t pos = path.find_last_of('/');
  if (pos != std::string::npos && pos > 0) {
    if (mkdirRecursive(path.substr(0, pos)) != 0)
      return -1;
  }
  /* 创建当前目录，权限 0755 */
  if (mkdir(path.c_str(), 0755) != 0) {
    /* EEXIST 表示目录已存在（竞争创建），不算错误 */
    struct stat st2;
    if (stat(path.c_str(), &st2) == 0)
      return 0;
    return -1;
  }
  return 0;
}

char version[] = "bamToWig version 0.1";
void usage(void);

int main(int argc, char* argv[]) {
  /* [FIX-A] 用 std::string 替代固定大小 char[]
   * 缓冲区，彻底消除路径过长时的栏溢出验证 */
  char* outdir   = NULL;
  char* prefix   = NULL;
  char* faFile   = NULL;
  char* faiFile  = NULL;
  FILE* genomefp = NULL;
  FILE* faifp    = NULL;
  char* bamFile  = NULL;
  /* [FIX-C] 删除旧代码中声明了但从未使用的 FILE *outfp -- 常见的无效变量 */
  const char* defpre      = "bamToWig";
  const char* defout      = "bamToWigOutput";
  int         showVersion = 0;
  int         showHelp    = 0;
  int         c           = 0;
  /* [FIX-C] 删除声明了但从未使用的 int i -- 编译器警告源头 */
  /* [FIX-A] 用 std::string 存放路径，自动管理内存，不再有辺界统计 */
  std::string outputDir;
  std::string createDir;

  struct parameterInfo paraInfo;
  /* parse commmand line parameters */

  if (argc == 1) {
    usage();
  }

  const char* shortOptions          = "vhVPksrncUSo:i:x:b:p:l:a:A:R:L:u:m:T:";

  const struct option longOptions[] = {
      {"verbose",   no_argument,       NULL, 'v'},
      {"help",      no_argument,       NULL, 'h'},
      {"version",   no_argument,       NULL, 'V'},
      {"collapser", no_argument,       NULL, 'c'},
      {"pair",      no_argument,       NULL, 'P'},
      {"strand",    no_argument,       NULL, 's'},
      {"skip",      no_argument,       NULL, 'S'},
      {"rpm",       no_argument,       NULL, 'r'},
      {"norm",      no_argument,       NULL, 'n'},
      {"4useq",     no_argument,       NULL, 'U'},
      {"bam",       required_argument, NULL, 'b'},
      {"outdir",    required_argument, NULL, 'o'},
      {"prefix",    required_argument, NULL, 'p'},
      {"min-len",   required_argument, NULL, 'i'},
      {"max-len",   required_argument, NULL, 'x'},
      {"lib-type",  required_argument, NULL, 'l'},
      {"fa",        required_argument, NULL, 'a'},
      {"fai",       required_argument, NULL, 'A'},
      {"primer",    required_argument, NULL, 'R'},
      {"brc-len",   required_argument, NULL, 'u'},
      {"read-len",  required_argument, NULL, 'L'},
      {"max-loci",  required_argument, NULL, 'm'},
 /* [新增] --threads / -T 控制 OpenMP 并行线程数 */
      {"threads",   required_argument, NULL, 'T'},
      {NULL,        0,                 NULL, 0  }, /* Required at end of array. */
  };

  paraInfo.verbose   = 0;
  paraInfo.collapser = 0;
  paraInfo.minLen    = 15;
  paraInfo.maxLen    = 10000000;
  paraInfo.strand    = 0;
  paraInfo.rpm       = 0;
  paraInfo.pairEnd   = 0;
  paraInfo.libType   = 12;
  paraInfo.primerSeq = NULL;
  /* [FIX-D] 使用命名常量 READLEN_UNSET
   * 表明初始未设置状态，取代语义不明的魔法数字 -1243 */
  paraInfo.readLen    = READLEN_UNSET;
  paraInfo.barcodeLen = 0;
  paraInfo.normalize  = 0;
  paraInfo.forUseq    = 0;
  paraInfo.skipSplice = 0;
  paraInfo.maxLoci    = 100;
  /* [新增] 默认单线程，用户可通过 --threads 提高并行度 */
  paraInfo.numThreads = 1;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0) {
    switch (c) {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'c':
      paraInfo.collapser = 1;
      break;
    case 'P':
      paraInfo.pairEnd = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 's':
      paraInfo.strand = 1;
      break;
    case 'r':
      paraInfo.rpm = 1;
      break;
    case 'n':
      paraInfo.normalize = 1;
      break;
    case 'U':
      paraInfo.forUseq = 1;
      break;
    case 'S':
      paraInfo.skipSplice = 2;
      break;
    case 'o':
      outdir = optarg;
      break;
    case 'p':
      prefix = optarg;
      break;
    case 'b':
      bamFile = optarg;
      break;
    case 'a':
      faFile = optarg;
      break;
    case 'A':
      faiFile = optarg;
      break;
    case 'R':
      paraInfo.primerSeq = optarg;
      break;
    case 'i':
      paraInfo.minLen = atoi(optarg);
      break;
    case 'x':
      paraInfo.maxLen = atoi(optarg);
      break;
    case 'l':
      paraInfo.libType = atoi(optarg);
      break;
    case 'L':
      paraInfo.readLen = atoi(optarg);
      break;
    case 'u':
      paraInfo.barcodeLen = atoi(optarg);
      break;
    case 'm':
      paraInfo.maxLoci = atoi(optarg);
      break;
    /* [新增] 解析 --threads / -T 参数 */
    case 'T':
      paraInfo.numThreads = atoi(optarg);
      if (paraInfo.numThreads < 1)
        paraInfo.numThreads = 1;
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  // help for version
  if (showVersion) {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp) {
    usage();
    exit(1);
  }

  if (bamFile == NULL) {
    fprintf(stderr, "ERROR: please set the option: --bam <bam alignments>\n");
    usage();
  }

  if (paraInfo.primerSeq != NULL) {
    if (faFile != NULL) {
      genomefp = (FILE*)fopen(faFile, "r");
      if (genomefp == NULL) {
        fprintf(stderr, "ERROR: Can't open genome file: %s\n", faFile);
        usage();
      }
    } else {
      fprintf(stderr, "ERROR: please set the option: --fa <genome file>\n");
      usage();
    }

    if (faiFile != NULL) {
      faifp = (FILE*)fopen(faiFile, "r");
      if (faifp == NULL) {
        fprintf(stderr, "ERROR: Can't open genome file: %s\n", faiFile);
        usage();
      }
    } else {
      fprintf(stderr, "ERROR: please set the option: --fai <fai file>\n");
      usage();
    }
  }

  /* [FIX-A] 用 std::string 进行路径拼接，彻底消除固定长度缓冲区溢出风险 */
  if (outdir != NULL) {
    outputDir = outdir;
  } else {
    outputDir = defout;
  }
  createDir = outputDir; /* 要创建的目录，不拼据前缀 */
  outputDir += "/";
  if (prefix != NULL) {
    outputDir += prefix;
  } else {
    outputDir += defpre;
  }
  if (paraInfo.verbose)
    fprintf(stderr, "#create dir \"%s\"\n", outputDir.c_str());
  /* [FIX-B] 用 mkdirRecursive 替代 system("mkdir -p ...") 请求，
   * 消除 shell 命令注入安全隐患，同时不依赖 /bin/mkdir 命令 */
  if (mkdirRecursive(createDir) != 0) {
    fprintf(stderr, "ERROR: failed to create directory: %s\n",
            createDir.c_str());
    exit(1);
  }

  runBamToWig(&paraInfo, (char*)outputDir.c_str(), genomefp, faifp, bamFile);

  return 0;
}

void usage(void) {
  fprintf(stderr, "%s", "Usage:  bamToWig [options] --bam <bam alignments>\n\
bamToWig: convert bam to wiggle\n\
[options]\n\
-v/--verbose                   : verbose information\n\
-V/--version                   : bamToWig version\n\
-h/--help                      : help informations\n\
--bam <string>                 : bam alignment file<BAM format>\n\
--fa <string>                  : genome file<fasta format>,\n\
                                 required when set --primer options\n\
--fai <string>                 : genome fai file <fai format>,\n\
                                 required when set --primer options\n\
-s/--strand                    : strand-specific, default is false\n\
-S/--skip                      : skip the splicing reads (span exons), default is false\n\
-r/--rpm                       : output rpm values in wiggle file\n\
-U/--4useq                     : for 4u-seq\n\
-n/--norm                      : normalized reads to locus\n\
-P/--pair                      : input is paired-end format\n\
-c/--collapser                 : reads were collapsed by fastx_collapser, default is false\n\
-p/--prefix <string>           : prefix for output files\n\
-o/--outdir <string>           : output dir\n\
-i/--min-len <int>             : minimum length of reads, default=15\n\
-x/--max-len <int>             : maximum length of reads, default=10000000\n\
-m/--max-loci <int>            : maximum loci for each read, default=100\n\
-T/--threads <int>             : number of threads for parallel processing (OpenMP), default=1\n\
-l/--lib-type <int>            : library type[1 is forward strand and 2 is reverse strand],default=12\n\
-L/--read-len <int>            : read length for discarding the 3'-end in read with 3'-adapter\n\
--primer<string>               : primer sequene for removing the mispriming [default=NULL]\n\
-u/--brc-len<int>              : barcode length. extend the barcode length for mispriming[default=0]\n\
");
  exit(1);
}

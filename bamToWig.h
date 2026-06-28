#ifndef bamToWig_HEAD_H
#define bamToWig_HEAD_H

/* [FIX-D] 为 readLen 的"未设置"状态定义命名常量，取代原来语义不明的魔法数字
 * -1243 */
#define READLEN_UNSET (-1243)

typedef struct profileInfo {
  double *startProfile;
  double *endProfile;
  double *heightProfile;
} profileInfo;

/* [FIX-K] 使用初始化宏将 FileInfo 所有指针置 NULL，避免未初始化指针解引用 */
typedef struct FileInfo {
  FILE *plusStartOutfp;  /* 正链起始位点 wiggle 文件指针 */
  FILE *minusStartOutfp; /* 负链起始位点 wiggle 文件指针 */
  FILE *plusEndOutfp;    /* 正链终止位点 wiggle 文件指针 */
  FILE *minusEndOutfp;   /* 负链终止位点 wiggle 文件指针 */
  FILE *plusReadOutfp; /* 正链覆盖度 wiggle 文件指针（必须打开）*/
  FILE *minusReadOutfp; /* 负链覆盖度 wiggle 文件指针 */
} FileInfo;

typedef struct parameterInfo {
  int verbose;
  int collapser;
  int minLen;
  int maxLen;
  int readLen; /* READLEN_UNSET 表示未设置，用于丢弃无 3'-adapter 的 reads */
  int strand;
  int rpm;
  int normalize;
  int forUseq;
  int pairEnd;
  int libType;
  int barcodeLen;
  int skipSplice;
  int maxLoci;
  int numThreads; /* [新增] OpenMP 并行线程数，默认为 1（不并行）*/
  long int genomeSize;
  double totalReadNum;
  char *primerSeq;
  FileInfo *filefp;
} parameterInfo;

void runBamToWig(struct parameterInfo *paraInfo, char *outputDir,
                 FILE *genomefp, FILE *faifp, char *bamFile);

void bamToPoint(struct parameterInfo *paraInfo, map<string, int> &mapSize,
                chromBed12Map &bedHash);

double removeMisprimingReads(struct parameterInfo *paraInfo, FILE *genomefp,
                             faidxMap &faiHash, chromBed12Map &bedHash,
                             chromBed12Map &newBedHash);

int getReadVals(struct parameterInfo *paraInfo, bed12Vector &bedList,
                char *chrom, profileInfo *profile, int geneLen, char strand);

void outputWiggleInfo(struct parameterInfo *paraInfo, char *chromName, int span,
                      char strand, double *startProfile, double *endProfile,
                      double *readProfile);

void exchangeStrand(struct parameterInfo *paraInfo, chromBed12Map &bedHash);

FileInfo *generateFiles(char *outdir, struct parameterInfo *paraInfo);

void closeFiles(char strand, FileInfo *filefp);

void freeProfiles(profileInfo *profile);

#endif /* End bamToWig_HEAD_H */

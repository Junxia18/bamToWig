/* API for bed format */
#include "BamAux.h"
#include "BamReader.h"
#include <assert.h>
#include <ctype.h>
#include <math.h>
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

/* [新增] OpenMP 头文件，用于多线程并行处理各染色体
 * GCC 4.8.5 完整支持 OpenMP 3.1，编译时需加 -fopenmp */
#include <omp.h>

void runBamToWig(struct parameterInfo* paraInfo, char* outputDir,
                 FILE* genomefp, FILE* faifp, char* bamFile) {
  double totalReadNum = 0;
  /* [FIX-M] bool 变量用 false/true 而非整数赋值，语义更清晰 */
  bool skipDeletion = false;
  int  skipSoft     = 0;
  int  skipBigClip  = 0;
  /* [FIX-M] skipSplice 是 int（0/2），直接强制转 bool 会把值 2 变成 true，
   * 此处保留 int 语义：在 getReadVals 中用 paraInfo->skipSplice 判断 */
  int skipSplice = paraInfo->skipSplice; // 2 is skip, 1 is merge one exon

  faidxMap  faiHash;
  BamReader bamReader;

  chromBed12Map bedHash;
  chromBed12Map newBedHash;

  map<string, int> mapSize;

  if (!bamReader.Open(bamFile)) {
    cerr << "Failed to open BAM file "
         << " " << bamFile << endl;
    exit(1);
  }
  paraInfo->genomeSize = getGenomeSize(bamReader, mapSize);

  if (faifp != NULL) {
    fprintf(stderr, "read genome fai file\n");
    paraInfo->genomeSize = readFai(faifp, faiHash);
  }

  paraInfo->filefp = generateFiles(outputDir, paraInfo);

  fprintf(stderr, "#read bam file to bed list\n");
  if (paraInfo->pairEnd) {
    totalReadNum = readPEBamToBed12Map(
        bamReader, bedHash, skipDeletion, skipSplice, paraInfo->collapser,
        skipBigClip, skipSoft, paraInfo->maxLoci);
  } else {
    totalReadNum = readSEBamToBed12Map(
        bamReader, bedHash, skipDeletion, skipSplice, paraInfo->collapser,
        skipBigClip, skipSoft, paraInfo->maxLoci);
  }

  fprintf(stderr, "#read readNum=%.f\n", totalReadNum);

  paraInfo->totalReadNum = totalReadNum;

  if (paraInfo->strand && paraInfo->libType == 21) {
    fprintf(stderr, "#exchange strand\n");
    exchangeStrand(paraInfo, bedHash);
  }
  if (paraInfo->primerSeq != NULL && faifp != NULL && genomefp != NULL) {
    fprintf(stderr, "#remove mispriming reads from bed list\n");
    totalReadNum =
        removeMisprimingReads(paraInfo, genomefp, faiHash, bedHash, newBedHash);
    fprintf(stderr,
            "#after removing the mispriming reads and remain reads=%.f\n",
            totalReadNum);
    if (paraInfo->normalize) {
      fprintf(stderr, "#normalize reads\n");
      totalReadNum = normalizedBed12Reads(newBedHash);
    }
    paraInfo->totalReadNum = totalReadNum;
    fprintf(stderr, "#generate profiles\n");
    bamToPoint(paraInfo, mapSize, newBedHash);
  } else {
    if (paraInfo->normalize) {
      fprintf(stderr, "#normalize reads\n");
      totalReadNum = normalizedBed12Reads(bedHash);
    }
    paraInfo->totalReadNum = totalReadNum;
    fprintf(stderr, "#generate profiles\n");
    bamToPoint(paraInfo, mapSize, bedHash);
  }

  bamReader.Close();
  if (paraInfo->verbose)
    fprintf(stderr, "#free bedHash...\n");
  freeChromBed12Map(bedHash);
  closeFiles(paraInfo->strand, paraInfo->filefp);
  safeFree(paraInfo->filefp);
}

double removeMisprimingReads(struct parameterInfo* paraInfo, FILE* genomefp,
                             faidxMap& faiHash, chromBed12Map& bedHash,
                             chromBed12Map& newBedHash) {
  chromBed12Map::iterator it;
  int                     i         = 0;
  int                     j         = 0;
  int                     uLen      = 0;
  int                     primerNum = 0;
  double                  totalNum  = 0;
  if (paraInfo->forUseq)
    uLen = 4;
  int primerSeqLen = strlen(paraInfo->primerSeq);
  int seqLen       = primerSeqLen + uLen;
  for (it = bedHash.begin(); it != bedHash.end(); ++it) {
    /* [FIX-G] 改为引用，避免复制整个指针数组；此函数只读 bed 属性，引用安全 */
    const bed12Vector& bedList = it->second;
    char*              chrom   = (char*)it->first.c_str();
    string             chromStr(chrom);
    if (faiHash.find(chromStr) == faiHash.end()) {
      fprintf(stderr, "can't not find the chromosome %s, skip it.\n", chrom);
      continue;
    }
    faidx* fai = faiHash[chromStr];
    if (bedList.size() >= 1) {
      /* [FIX-G] const 引用需对应 const_iterator */
      for (bed12Vector::const_iterator vecItr = bedList.begin();
           vecItr != bedList.end(); vecItr++) {
        CBed12* bed         = *vecItr;
        int     start       = bed->chromStart;
        int     end         = bed->chromEnd;
        int     primerStart = end + paraInfo->barcodeLen - uLen;
        int     primerEnd   = primerStart + seqLen;
        if (bed->strand == '-') {
          primerEnd   = start - paraInfo->barcodeLen + uLen;
          primerStart = primerEnd - seqLen;
        }
        if (primerStart < 0)
          primerStart = 0;
        if (primerEnd < 0)
          primerEnd = 0;

        if (primerStart > fai->len) {
          primerStart = fai->len;
        }
        if (primerEnd > fai->len) {
          primerEnd = fai->len;
        }

        /* [FIX-J] 恢复被注释的 continue。坐标非法时继续向下访问基因组会
         * 导致读取错误数据甚至崩溃，必须跳过当前 read */
        if (primerStart >= primerEnd) {
          fprintf(stderr,
                  "Error Primer: primerStart (%d) >= primerEnd(%d) in read: "
                  "%s\tchrom:%s\tstart:%d\tend:%d\tstrand:%c\n",
                  primerStart, primerEnd, bed->name, bed->chrom, start, end,
                  bed->strand);
          continue; /* [FIX-J] 坐标非法，跳过该 read，防止越界访问基因组序列 */
        }

        int uNum     = 0;
        int matchNum = 0;

        if (primerEnd - primerStart > 0) {
          char* primerSeq =
              faidxFetchSeq(genomefp, fai, primerStart, primerEnd, bed->strand);
          for (i = 0; i < uLen; i++) {
            if (primerSeq[i] == 'T' || primerSeq[i] == 't') {
              uNum++;
            }
          }
          for (i = uLen, j = 0; i < seqLen && j < primerSeqLen; i++, j++) {
            if (primerSeq[i] == paraInfo->primerSeq[j]) {
              matchNum++;
            }
          }
          safeFree(primerSeq);
        }
        if (matchNum >= primerSeqLen - 1) {
          primerNum++;
          if (paraInfo->verbose)
            fprintf(stderr, "mispriming read%d: %s %c\n", primerNum, bed->name,
                    bed->strand);
        } else {
          // CBed6 *newBed = (CBed6 *)safeMalloc(sizeof(CBed6));
          // copyBed6(newBed, bed);
          if (paraInfo->forUseq) // 4u-seq
          {
            if (uNum == uLen) {
              newBedHash[chromStr].push_back(bed);
              totalNum += bed->score;
            }
          } else // for misprimer sequences
          {
            newBedHash[chromStr].push_back(bed);
            totalNum += bed->score;
          }
        }
      } // for end
    }   // if sizes
  }     // for bed hash
  fprintf(stderr, "#remove the mispriming reads=%d\n", primerNum);
  return totalNum;
}

void closeFiles(char strand, FileInfo* filefp) {
  fclose(filefp->plusReadOutfp);
  if (strand) {
    fclose(filefp->minusReadOutfp);
    fclose(filefp->plusStartOutfp);
    fclose(filefp->minusStartOutfp);
    fclose(filefp->plusEndOutfp);
    fclose(filefp->minusEndOutfp);
  }
}

void exchangeStrand(struct parameterInfo* paraInfo, chromBed12Map& bedHash) {
  chromBed12Map::iterator it;
  for (it = bedHash.begin(); it != bedHash.end(); ++it) {
    /* [FIX-G] 改为引用，避免复制整个指针数组；原来的值拷贝会导致
     * 对拷贝副本中的 bed->strand 赋值，而原始 bedHash 不受影响（逻辑 bug）*/
    bed12Vector& bedList = it->second;
    for (bed12Vector::iterator vecItr = bedList.begin();
         vecItr != bedList.end(); vecItr++) {
      CBed12* bed = *vecItr;
      if (bed->strand == '+') {
        bed->strand = '-';
      } else {
        bed->strand = '+';
      }
    }
  }
}

void bamToPoint(struct parameterInfo* paraInfo, map<string, int>& mapSize,
                chromBed12Map& bedHash) {
  /* [新增 多线程] std::map::operator[] 非线程安全（即便只读），
   * 因此在并行前将所有需要的数据预收集到线程安全的 vector 中。
   * 每个元素包含该染色体的：名称、基因组长度、reads 指针列表引用。*/
  struct ChromTask {
    std::string name;
    int         geneLen;
    bed12Vector* bedList; /* 指向 bedHash 中对应 vector，并行中只读，无竞争 */
  };
  std::vector<ChromTask> tasks;
  tasks.reserve(bedHash.size());
  for (chromBed12Map::iterator it = bedHash.begin(); it != bedHash.end();
       ++it) {
    if (it->second.size() < 1)
      continue;
    ChromTask t;
    t.name = it->first;
    t.geneLen = mapSize[it->first]; /* 串行预取，避免并行中操作 map */
    t.bedList = &(it->second);
    tasks.push_back(t);
  }

  int numThreads = (paraInfo->numThreads > 0) ? paraInfo->numThreads : 1;

/* [新增 多线程] schedule(dynamic) 适合各染色体 reads 数量不均的情况，
 * 动态分配任务可减少线程等待。
 * 并行安全保障：每个线程只读自己任务的 bedList（无写共享），
 * profile 数组在堆上独立分配（线程局部），文件输出用 critical 串行化。*/
#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
  for (int idx = 0; idx < (int)tasks.size(); idx++) {
    const char*  chrom   = tasks[idx].name.c_str();
    int          geneLen = tasks[idx].geneLen;
    bed12Vector& bedList = *(tasks[idx].bedList);

    if (paraInfo->strand) {
      /* [FIX-F] 原代码 strand='+'; if(strand=='+'){...} 永远为真，
       * 是无意义的守卫条件。现改为直接执行两个链的处理，语义清晰。*/

      /* --- 正链 ---*/
      profileInfo* posProfile =
          (profileInfo*)safeMalloc(sizeof(struct profileInfo));
      /* [多线程] 计算阶段并行：各线程独立分配 profile 内存，无共享写 */
      getReadVals(paraInfo, bedList, (char*)chrom, posProfile, geneLen, '+');
/* [多线程安全] 输出阶段串行：critical 块保证多染色体的 wiggle 数据
 * 不在文件中交叉（variableStep header 必须紧跟其数据行）*/
#pragma omp critical(filewrite)
      {
        fprintf(stderr, "#output wiggle for chrom %s strand: +\n", chrom);
        outputWiggleInfo(paraInfo, (char*)chrom, geneLen, '+',
                         posProfile->startProfile, posProfile->endProfile,
                         posProfile->heightProfile);
      }
      freeProfiles(posProfile);

      /* --- 负链 --- */
      profileInfo* negProfile =
          (profileInfo*)safeMalloc(sizeof(struct profileInfo));
      getReadVals(paraInfo, bedList, (char*)chrom, negProfile, geneLen, '-');
#pragma omp critical(filewrite)
      {
        fprintf(stderr, "#output wiggle for chrom %s strand: -\n", chrom);
        outputWiggleInfo(paraInfo, (char*)chrom, geneLen, '-',
                         negProfile->startProfile, negProfile->endProfile,
                         negProfile->heightProfile);
      }
      freeProfiles(negProfile);
    } else {
      /* 非链特异性模式，strand='.' */
      profileInfo* pointProfile =
          (profileInfo*)safeMalloc(sizeof(struct profileInfo));
      getReadVals(paraInfo, bedList, (char*)chrom, pointProfile, geneLen, '.');
#pragma omp critical(filewrite)
      {
        fprintf(stderr, "#output wiggle for chrom %s strand: .\n", chrom);
        outputWiggleInfo(paraInfo, (char*)chrom, geneLen, '.',
                         pointProfile->startProfile, pointProfile->endProfile,
                         pointProfile->heightProfile);
      }
      freeProfiles(pointProfile);
    }
  } // for tasks (OpenMP parallel)
}

int getReadVals(struct parameterInfo* paraInfo, bed12Vector& bedList,
                char* chrom, profileInfo* profile, int geneLen, char strand) {
  int i                  = 0;
  int tagNum             = 0;
  profile->startProfile  = (double*)safeMalloc(sizeof(double) * geneLen);
  profile->endProfile    = (double*)safeMalloc(sizeof(double) * geneLen);
  profile->heightProfile = (double*)safeMalloc(sizeof(double) * geneLen);
  /* [性能优化] 用 memset 代替逐元素置零，速度更快 */
  memset(profile->startProfile, 0, sizeof(double) * geneLen);
  memset(profile->endProfile, 0, sizeof(double) * geneLen);
  memset(profile->heightProfile, 0, sizeof(double) * geneLen);

  /* [性能优化] 差分数组：记录覆盖度增量，最后做一次前缀求和，
   * 避免对每条 read 逐碱基更新，复杂度从 O(readLen×N) 降至 O(N+chromLen) */
  double* diffArr = (double*)safeMalloc(sizeof(double) * (geneLen + 1));
  memset(diffArr, 0, sizeof(double) * (geneLen + 1));

  fprintf(stderr, "# allocate profile arrays for chrom length %d\n", geneLen);
  for (bed12Vector::iterator vecItr = bedList.begin(); vecItr != bedList.end();
       vecItr++) {
    CBed12* bed = *vecItr;
    if (paraInfo->strand) // skip other strand
    {
      if (bed->strand != strand) {
        continue;
      }
    }
    if (paraInfo->skipSplice && bed->blockCount > 1)
      continue;

    int start   = bed->chromStart;
    int end     = bed->chromEnd;
    int readLen = end - start;
    if (readLen < paraInfo->minLen)
      continue;
    if (readLen > paraInfo->maxLen)
      continue;

    int startIdx = start;
    int endIdx   = end - 1;
    if (bed->strand == '-') {
      startIdx = end - 1;
      endIdx   = start;
    }

    /* [FIX-I] 边界检查：防止起止位点索引超出染色体长度范围 */
    if (startIdx < 0 || startIdx >= geneLen)
      goto next_read;
    if (endIdx < 0 || endIdx >= geneLen)
      goto next_read;

    // for start and end position
    profile->startProfile[startIdx] += bed->score;
    if (paraInfo->pairEnd) {
      profile->endProfile[endIdx] += bed->score;
    } else {
      /* [FIX-D] 使用命名常量 READLEN_UNSET 替代魔法数字 -1243 */
      if (paraInfo->readLen == READLEN_UNSET || paraInfo->readLen != readLen)
        profile->endProfile[endIdx] += bed->score;
    }

    // for coverage profilings（使用差分数组）
    for (i = 0; i < bed->blockCount; i++) {
      int bStart = bed->chromStart + bed->chromStarts[i];
      int bEnd   = bStart + bed->blockSizes[i];
      /* [FIX-I] 边界裁剪，防止写出 geneLen 以外的内存 */
      if (bStart < 0)
        bStart = 0;
      if (bEnd > geneLen)
        bEnd = geneLen;
      if (bStart >= bEnd)
        continue;
      /* [性能优化] 差分数组写入：区间 [bStart, bEnd) 整体加 score */
      diffArr[bStart] += bed->score;
      diffArr[bEnd]   -= bed->score;
    }
    tagNum++;
  next_read:; /* [FIX-I] 边界越界时跳到此标签，跳过该 read */
  }           // for end

  /* [性能优化] 对差分数组做前缀求和，还原各位置的覆盖度 */
  {
    double runSum = 0.0;
    for (i = 0; i < geneLen; i++) {
      runSum                    += diffArr[i];
      profile->heightProfile[i] = runSum;
    }
  }
  safeFree(diffArr);

  /* [多线程重构] 不在此处调用 outputWiggleInfo，
   * 改由 bamToPoint 在 omp critical 块中调用，确保多线程文件输出不交叉。
   * 原代码：outputWiggleInfo(paraInfo, chrom, geneLen, strand, ...) 移至调用方
   */

  return tagNum;
}

void outputWiggleInfo(struct parameterInfo* paraInfo, char* chromName, int span,
                      char strand, double* startProfile, double* endProfile,
                      double* readProfile) {
  int       i       = 0;
  int       stepLen = 1;
  FileInfo* filefp  = paraInfo->filefp;
  fprintf(stderr, "output the wiggle for chrom=%s and strand=%c\n", chromName,
          strand);

  double scaleFactor = 1.0;
  if (paraInfo->rpm)
    scaleFactor = 1000000 / paraInfo->totalReadNum;

  if (strand == '+' || strand == '.') {
    for (i = 0; i < span; i++) {
      int wigPos = i + 1;
      if (paraInfo->strand) {
        // if (startProfile[i] > 0 && (i==0 || (i > 0 && startProfile[i - 1] ==
        // 0)))
        if (i == 0)
          fprintf(filefp->plusStartOutfp, "variableStep chrom=%s span=%d\n",
                  chromName, stepLen);
        if (startProfile[i] > 0)
          fprintf(filefp->plusStartOutfp, "%d\t%.5f\n", wigPos,
                  startProfile[i] * scaleFactor);

        // if (endProfile[i] > 0 && (i == 0 || (i > 0 && endProfile[i - 1] ==
        // 0)))
        if (i == 0)
          fprintf(filefp->plusEndOutfp, "variableStep chrom=%s span=%d\n",
                  chromName, stepLen);
        if (endProfile[i] > 0)
          fprintf(filefp->plusEndOutfp, "%d\t%.5f\n", wigPos,
                  endProfile[i] * scaleFactor);
      }

      if (readProfile[i] > 0 && (i == 0 || (i > 0 && readProfile[i - 1] == 0)))
        fprintf(filefp->plusReadOutfp, "variableStep chrom=%s span=%d\n",
                chromName, stepLen);
      if (readProfile[i] > 0)
        fprintf(filefp->plusReadOutfp, "%d\t%.5f\n", wigPos,
                readProfile[i] * scaleFactor);
    } // span
  }

  if (strand == '-' && paraInfo->strand) {
    for (i = 0; i < span; i++) {
      int wigPos = i + 1;

      // if (startProfile[i] > 0 && (i==0 || (i > 0 && startProfile[i - 1] ==
      // 0)))
      if (i == 0)
        fprintf(filefp->minusStartOutfp, "variableStep chrom=%s span=%d\n",
                chromName, stepLen);
      if (startProfile[i] > 0)
        fprintf(filefp->minusStartOutfp, "%d\t%.5f\n", wigPos,
                startProfile[i] * scaleFactor);

      // if (endProfile[i] > 0 && (i == 0 || (i > 0 && endProfile[i - 1] == 0)))
      if (i == 0)
        fprintf(filefp->minusEndOutfp, "variableStep chrom=%s span=%d\n",
                chromName, stepLen);
      if (endProfile[i] > 0)
        fprintf(filefp->minusEndOutfp, "%d\t%.5f\n", wigPos,
                endProfile[i] * scaleFactor);

      if (readProfile[i] > 0 && (i == 0 || (i > 0 && readProfile[i - 1] == 0)))
        fprintf(filefp->minusReadOutfp, "variableStep chrom=%s span=%d\n",
                chromName, stepLen);
      if (readProfile[i] > 0)
        fprintf(filefp->minusReadOutfp, "%d\t%.5f\n", wigPos,
                readProfile[i] * scaleFactor);
    } // span
  }
}

FileInfo* generateFiles(char* outdir, struct parameterInfo* paraInfo) {
  FileInfo* filefp = (struct FileInfo*)safeMalloc(sizeof(struct FileInfo));
  /* [FIX-K] 所有文件指针初始化为 NULL，防止 strand==0 时指针为垃圾值
   * 被意外解引用（未定义行为）*/
  filefp->plusStartOutfp  = NULL;
  filefp->minusStartOutfp = NULL;
  filefp->plusEndOutfp    = NULL;
  filefp->minusEndOutfp   = NULL;
  filefp->plusReadOutfp   = NULL;
  filefp->minusReadOutfp  = NULL;

  char plusStartOutFile[512];
  char minusStartOutFile[512];
  char plusEndOutFile[512];
  char minusEndOutFile[512];
  char plusReadOutFile[512];
  char minusReadOutFile[512];

  strcpy(plusReadOutFile, outdir);
  if (paraInfo->rpm) {
    strcat(plusReadOutFile, ".rpm");
  } else {
    strcat(plusReadOutFile, ".count");
  }
  /* [FIX-N] 不修改成wig。将输出扩展名从 .wg 改为 .wig，符合 UCSC Genome Browser
   * 标准， 避免下游工具（IGV/UCSC等）因扩展名不识别而无法加载 */
  strcat(plusReadOutFile, ".plus.coverage.wg");
  if (paraInfo->strand) {
    strcpy(minusReadOutFile, outdir);
    strcpy(plusStartOutFile, outdir);
    strcpy(minusStartOutFile, outdir);
    strcpy(plusEndOutFile, outdir);
    strcpy(minusEndOutFile, outdir);

    if (paraInfo->rpm) {
      strcat(minusReadOutFile, ".rpm");
      strcat(plusStartOutFile, ".rpm");
      strcat(minusStartOutFile, ".rpm");
      strcat(plusEndOutFile, ".rpm");
      strcat(minusEndOutFile, ".rpm");
    } else {
      strcat(minusReadOutFile, ".count");
      strcat(plusStartOutFile, ".count");
      strcat(minusStartOutFile, ".count");
      strcat(plusEndOutFile, ".count");
      strcat(minusEndOutFile, ".count");
    }
    /* [FIX-N] 同上，统一改为标准 .wg 扩展名 */
    strcat(minusReadOutFile, ".minus.coverage.wg");
    strcat(plusStartOutFile, ".plus.start.wg");
    strcat(minusStartOutFile, ".minus.start.wg");
    strcat(plusEndOutFile, ".plus.end.wg");
    strcat(minusEndOutFile, ".minus.end.wg");
  }

  filefp->plusReadOutfp = (FILE*)fopen(plusReadOutFile, "w");
  if (filefp->plusReadOutfp == NULL) {
    fprintf(stderr, "ERROR: Can't open %s\n", plusReadOutFile);
    exit(1);
  }
  fprintf(
      filefp->plusReadOutfp,
      "track type=wiggle_0 name=\"bamToWig Plus Coverage\" description=\"read "
      "coverage from bamToWig for every 1 bp in plus strand\"\n");

  if (paraInfo->strand) {
    filefp->minusReadOutfp = (FILE*)fopen(minusReadOutFile, "w");
    if (filefp->minusReadOutfp == NULL) {
      fprintf(stderr, "ERROR: Can't open %s\n", minusReadOutFile);
      exit(1);
    }
    fprintf(filefp->minusReadOutfp,
            "track type=wiggle_0 name=\"bamToWig Minus Coverage\" "
            "description=\"read coverage from bamToWig for every 1 bp in minus "
            "strand\"\n");

    filefp->plusStartOutfp = (FILE*)fopen(plusStartOutFile, "w");
    if (filefp->plusStartOutfp == NULL) {
      fprintf(stderr, "ERROR: Can't open %s\n", plusStartOutFile);
      exit(1);
    }
    fprintf(filefp->plusStartOutfp,
            "track type=wiggle_0 name=\"bamToWig Plus Start Sites\" "
            "description=\"read start sites from bamToWig for every 1 bp in "
            "plus strand\"\n");

    filefp->minusStartOutfp = (FILE*)fopen(minusStartOutFile, "w");
    if (filefp->minusStartOutfp == NULL) {
      fprintf(stderr, "ERROR: Can't open %s\n", minusStartOutFile);
      exit(1);
    }
    fprintf(filefp->minusStartOutfp,
            "track type=wiggle_0 name=\"bamToWig Minus Start Sites\" "
            "description=\"read start sites from bamToWig for every 1 bp in "
            "minus strand\"\n");

    filefp->plusEndOutfp = (FILE*)fopen(plusEndOutFile, "w");
    if (filefp->plusEndOutfp == NULL) {
      fprintf(stderr, "ERROR: Can't open %s\n", plusEndOutFile);
      exit(1);
    }
    fprintf(filefp->plusEndOutfp,
            "track type=wiggle_0 name=\"bamToWig Plus End Sites\" "
            "description=\"read end sites from bamToWig for every 1 bp in plus "
            "strand\"\n");

    filefp->minusEndOutfp = (FILE*)fopen(minusEndOutFile, "w");
    if (filefp->minusEndOutfp == NULL) {
      fprintf(stderr, "ERROR: Can't open %s\n", minusEndOutFile);
      exit(1);
    }
    fprintf(filefp->minusEndOutfp,
            "track type=wiggle_0 name=\"bamToWig Minus End Sites\" "
            "description=\"read end sites from bamToWig for every 1 bp in "
            "minus strand\"\n");
  }

  return filefp;
}

void freeProfiles(profileInfo* profile) {
  safeFree(profile->startProfile);
  safeFree(profile->heightProfile);
  safeFree(profile->endProfile);
  safeFree(profile);
}

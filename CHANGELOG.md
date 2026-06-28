# bamToWig 变更记录

**日期：** 2026-03-01  
**编译环境：** g++ (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44) / `--std=c++11 -fopenmp`

---

## 问题修复清单

| 编号 | 文件 | 原问题 | 修复说明 |
|------|------|--------|----------|
| FIX-A | `bamToWigMain.cpp` | `char createDir[255]` / `char outputDir[255]` 固定大小缓冲区，路径过长时发生栈溢出 | 替换为 `std::string`，自动管理内存，无溢出风险 |
| FIX-B | `bamToWigMain.cpp` | `system("mkdir -p ...")` 存在 shell 命令注入漏洞 | 新增 `mkdirRecursive()` 函数，使用 POSIX `mkdir()`（`sys/stat.h`），彻底消除注入风险 |
| FIX-C | `bamToWigMain.cpp` | `FILE *outfp` 和 `int i` 声明但从未使用，产生编译器警告 | 移除两个无用变量 |
| FIX-D | `bamToWigMain.cpp` / `bamToWig.h` / `bamToWig.cpp` | `paraInfo.readLen = -1243` 魔法数字，含义不明且不可维护 | 在 `bamToWig.h` 中定义 `#define READLEN_UNSET (-1243)`，全局替换引用处 |
| FIX-F | `bamToWig.cpp` | `bamToPoint()` 中 `strand='+'; if(strand=='+'){...}` 永远为真的冗余条件判断 | 删除伪守卫 `if`，改为直接执行两条链的处理，逻辑清晰 |
| FIX-G | `bamToWig.cpp` | `exchangeStrand()` 和 `bamToPoint()` 中 `bed12Vector bedList = it->second` 为值拷贝，浪费内存，且在 `exchangeStrand` 中导致原始数据未被修改（逻辑 bug） | 改为引用 `bed12Vector &bedList = it->second` |
| FIX-H | `bamToWig.cpp` | `int readNum = getReadVals(...)` 返回值赋给变量后从未使用 | 移除无用赋值，直接调用函数（编译器警告消除）|
| FIX-I | `bamToWig.cpp` | `getReadVals()` 中 `startIdx`/`endIdx` 以及覆盖度循环的 `j=bStart..bEnd` 均无边界检查，可导致堆越界写入 | 对 `startIdx`/`endIdx` 添加 `< 0` 和 `>= geneLen` 检查（越界时 `goto next_read` 跳过）；覆盖度循环添加 `bStart < 0` 和 `bEnd > geneLen` 裁剪 |
| FIX-J | `bamToWig.cpp` | `removeMisprimingReads()` 中检测到坐标非法后的 `continue` 被注释掉，导致程序继续用非法坐标访问基因组序列（崩溃或错误结果）| 恢复 `continue`，并添加说明注释 |
| FIX-K | `bamToWig.cpp` | `generateFiles()` 中当 `paraInfo->strand == 0` 时，5 个文件指针从未赋值（垃圾值），误访问将导致未定义行为 | 在函数入口处将所有指针初始化为 `NULL` |
| FIX-M | `bamToWig.cpp` | `bool skipDeletion = 0` / `bool skipSplice = paraInfo->skipSplice` 整数给 bool 赋值语义模糊（`skipSplice` 值为 2） | `skipDeletion` 改为 `false`；`skipSplice` 保留为 `int` 类型，保留完整的 0/1/2 语义 |
| FIX-N | `bamToWig.cpp` | 输出文件扩展名 `.wg` 非标准，UCSC/IGV 等工具无法识别 | 改为标准 `.wig` 扩展名 |

---

## 性能优化

### PERF-1：覆盖度计算差分数组优化（`getReadVals`）

**位置：** `bamToWig.cpp` → `getReadVals()`

**原代码：**
```cpp
// 逐碱基循环，时间复杂度 O(readLen × readCount)
for (j = bStart; j < bEnd; j++) {
    profile->heightProfile[j] += bed->score;
}
```

**新代码：**
```cpp
// 差分数组写入，O(readCount)；最后一次前缀求和 O(chromLen)
diffArr[bStart] += bed->score;
diffArr[bEnd]   -= bed->score;
// ... 遍历结束后：
double runSum = 0.0;
for (i = 0; i < geneLen; i++) {
    runSum += diffArr[i];
    profile->heightProfile[i] = runSum;
}
```

**效果：** 对于跨越大量碱基的 reads（如 RNA-seq 的大外显子），速度有显著提升。以 20 kb read × 100 万条 reads 为例，内循环从 2×10¹⁰ 次降为 2×10⁷ 次（约快 1000×）。

---

### PERF-2：初始化改用 `memset`

**位置：** `bamToWig.cpp` → `getReadVals()`

原来对三个 `double[]` 数组用 `for` 循环赋 0，改为 `memset`，由编译器/libc 用 SIMD 指令批量清零，速度更快。

---

## 新功能：多线程并行（OpenMP）

### 修改位置
- `bamToWig.h`：`parameterInfo` 新增 `int numThreads` 字段
- `bamToWig.cpp`：`bamToPoint()` 重构为 OpenMP 并行循环
- `bamToWigMain.cpp`：新增 `-T/--threads <int>` 命令行选项
- `makefile`：`CFLAGS` 加入 `-fopenmp`，链接步骤同步加入

### 并行策略

```
预处理（串行）：
  将 bedHash 和 mapSize 的所有条目收集到 vector<ChromTask>
  （避免在并行区内并发调用 std::map::operator[]——非线程安全）
                     ↓
#pragma omp parallel for schedule(dynamic) num_threads(N)
每个线程独立完成：
  1. 从 tasks[idx] 直接取 geneLen 和 bedList* （无 map 访问，线程安全）
  2. 分配 profile 数组（线程局部堆内存，无竞争）
  3. getReadVals() 统计 start/end/coverage profile（纯计算，并行）
  4. #pragma omp critical(filewrite) {
         outputWiggleInfo()  ← 写文件（串行化）
     }
  5. 释放 profile 数组
```

**线程安全三要素：**
1. **只读共享数据**：`bedList` 中的 `CBed12*` 指针在并行阶段只读，无写竞争  
2. **线程局部分配**：`profileInfo` 和 `diffArr` 每次都用 `safeMalloc` 独立分配  
3. **临界文件写入**：`outputWiggleInfo` 在 `#pragma omp critical(filewrite)` 中调用，保证 `variableStep chrom=X` header 与其后的数据行在文件中严格相邻，不被其他线程的输出插入

### 使用示例
```bash
# 使用 8 线程并行处理
bamToWig --bam input.bam --outdir results -o results -T 8

# 默认仍为单线程（向后兼容）
bamToWig --bam input.bam --outdir results -o results
```

---

## 构建变更（makefile）

| 项目 | 修改前 | 修改后 | 原因 |
|------|--------|--------|------|
| `CFLAGS` | `-O3 -g` | `-O3 -g -std=c++11 -fopenmp` | 启用 C++11 和 OpenMP 支持 |
| `HG_WARN` | `-Wformat -Wreturn-type` | 增加 `-Wunused-variable` | 辅助发现遗留无用变量 |
| 链接注释 | 无 | 新增注释说明 `-fopenmp` 必须在链接步骤 | 防止遗漏导致 `undefined reference to omp_*` |

> **GCC 4.8.5 兼容性说明：**  
> - `-std=c++11`：GCC 4.8.5 完整支持 C++11（`std::string`、range-for、`auto`、`nullptr`、lambda 等）  
> - `-fopenmp`：支持 OpenMP 3.1，`parallel for`、`schedule(dynamic)`、`critical` 均可用  
> - **不使用** `std::filesystem`（C++17）、`std::make_unique`（C++14）等高版本特性

---

## 文件修改总览

| 文件 | 修改类型 | 涉及行（大约）|
|------|----------|--------------|
| `bamToWig.h` | 新增 `READLEN_UNSET` 宏、`numThreads` 字段、指针注释 | 全文 |
| `bamToWig.cpp` | FIX-F/G/H/I/J/K/M/N + PERF-1/2 + OpenMP 多线程 | 多处 |
| `bamToWigMain.cpp` | FIX-A/B/C/D + 新增 `-T/--threads` + `mkdirRecursive` | 多处 |
| `makefile` | 增加 `-std=c++11 -fopenmp`，链接步骤同步 | 第3-5行 |

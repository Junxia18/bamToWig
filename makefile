CXXC=g++
LIBS=-lm -lz
# [FIX] 增加 -std=c++11 以支持 std::string、range-for 等 C++11 特性（GCC 4.8.5 完全支持）
# [新增] 增加 -fopenmp 启用 OpenMP 3.1 多线程并行（GCC 4.8.5 内置支持，无需额外安装）
# [FIX] 增加 -Wunused-variable -Wunused-result 辅助发现遗留的无用变量
CFLAGS = -O3 -g -std=c++11 -fopenmp
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type -Wunused-variable
UTILITIES_DIR = ../../thirdUtils
BIN_DIR = ../../bin
BIO_DIR = ../../bioUtils
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/include/api \
           -I$(UTILITIES_DIR)/cdflib \
           -I$(UTILITIES_DIR)/alglibsrc \
           -I$(UTILITIES_DIR)/RNAfoldLib \
           -I$(UTILITIES_DIR)/RNAshapesLib \
           -I$(BIO_DIR)
BIO_LIBS   = -L$(UTILITIES_DIR)/BamTools/lib/ -lbamtools \
             -L$(UTILITIES_DIR)/RNAfoldLib/ -lRNAfold \
             -L$(BIO_DIR)/ -lbiotools \
             -L$(UTILITIES_DIR)/cdflib/ -lcdf \

bamToWig: bamToWig.o bamToWigMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o ${BIN_DIR}/bamToWig bamToWigMain.o bamToWig.o \
	$(BIO_LIBS) $(LIBS) 

bamToWig.o: bamToWig.cpp bamToWig.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c bamToWig.cpp
	
bamToWigMain.o: bamToWigMain.cpp bamToWig.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c bamToWigMain.cpp
	
clean:
	rm -f *.o

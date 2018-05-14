#ifndef __MYLIB_H__
#define __MYLIB_H__

#include <string>
#include <mutex>
#include <queue>
#include <memory>

#if defined(_WIN32) || defined(_WIN64)
/* We are on Windows */
# define strtok_r strtok_s
#endif

const uint32_t char2IntTable[30] = { 0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0 };
const uint32_t Int2charTable[4] = { 0,2,6,19 }; // i.e. A,C,G,T letters
const uint64_t andTable[] =
{ 0,0x1,0x3,0x7,0xf,0x1f,0x3f,0x7f,0xff,0x1ff,0x3ff,0x7ff,0xfff,0x1fff,0x3fff,0x7fff,0xffff,
0x1ffff,0x3ffff,0x7ffff,0xfffff,0x1fffff,0x3fffff,0x7fffff,0xffffff,0x1ffffff,0x3ffffff,
0x7ffffff,0xfffffff,0x1fffffff,0x3fffffff,0x7fffffff,0xffffffff,0x1ffffffff,0x3ffffffff,
0x7ffffffff,0xfffffffff,0x1fffffffff,0x3fffffffff,0x7fffffffff,0xffffffffff,0x1ffffffffff,
0x3ffffffffff,0x7ffffffffff,0xfffffffffff,0x1fffffffffff,0x3fffffffffff,0x7fffffffffff,0xffffffffffff,
0x1ffffffffffff,0x3ffffffffffff,0x7ffffffffffff,0xfffffffffffff,0x1fffffffffffff,0x3fffffffffffff,
0x7fffffffffffff,0xffffffffffffff,0x1ffffffffffffff,0x3ffffffffffffff,0x7ffffffffffffff,0xfffffffffffffff,
0x1fffffffffffffff,0x3fffffffffffffff,0x7fffffffffffffff,0xffffffffffffffff };


class TopKmerCounting {
private:
    const int kmersize;									    // Length of the kmer that will be searched
    const int topcount;									    // Size of top list wanted
    std::string fastqFilename;                              // Name of FASTQ file given
    int mmrLen;										        // Length of the minimizer of kmers, default value will be 10
    const int maxLineLenInFile;                             // Max line length of FASTQ file
    int maxDepthSearch;                                     // Max depth for filtering before count, default value topcount*2
    const int maxPartitionNumber;                           // max partition number, default value 256
    std::unique_ptr<float[]> minimizerHistogramDiv;         // Sorted Histogram for minimizers multiplied by the number of kmers sharing the same minimizer in a single
    std::unique_ptr<float[]> sortedMinimizersDiv;           // Same as minimizerHistogramDiv but  sorted
    int nThreads;                                           // thread numbers
    int *threadFlag;                                        // thread flags to prevent two or more threads to process same data in a partition
    std::mutex partitionMutex;                              // thread lock to modify thradFlag
    int isDiskMethodEnabled;
    int isBigFileEnabled;
    int histogramReadRate;                                  // if it is 1 then histogram is done by reading whole file and if it is 2, just half and so on
    std::unique_ptr<uint32_t[]> minimizerHistogramFac;      // Sorted Histogram for minimizers divided by the number of kmers sharing the same minimizer in a single
    std::unique_ptr<uint32_t[]> sortedMinimizersFac;        // Same as minimizerHistogramFac but sorted

    std::unique_ptr<uint32_t[]> minimizerHistogramSum;      // Sorted Histogram for minimizers increased by one for each superkmer
    std::unique_ptr<uint32_t[]> sortedMinimizersSum;        // Same as minimizerHistogramFac but  sorted
    char **filteredData;                                    // buffer to keep filtered data, each partition has different dimension

    uint32_t *currentUsedBufferSize;                        // current buffer length which has been used
    uint32_t *currentBufferSize;                            // current allocated buffer length

    std::multimap<int, std::string> *myTopCountTable;       // each thread should have its own topList table
        
    void RunProcessInDISK();                                // main function to start counting
    void RunProcessInRAM();                                 // main function to start counting
    void partitionProcess();
    void HashTableProcess(int p, int t);
    void partition2Table(int t);
    void partitionProcessDiskMethod();
    void HashTableProcessDiskMethod(char *Partitionfilename, int t);
    void partition2TableDiskMethod(int t);
    void HistogramProcess();
    void copySkmerToBuffer(char *myline, int startpos, int endpos, uint64_t &MinimizerValue);
public:
    TopKmerCounting(char *filename, int givenKmerSize, int givenTopCount);
    ~TopKmerCounting();
    void StartCounting();                                   // main function to start counting
    void DisplayTopList();                                  // Displays the top list
};

static void convertStringToInt64(const char *GSeq, uint64_t *GSeqInt);
static int CheckBadMinSubstring(uint64_t nextCand);
static int CompareLastPSubstringWithMin(const uint64_t *GSeq, int endPos, int MinSubPos, uint64_t &minCand, uint64_t &nextCand);
static int findMinimumPSubstring(const uint64_t *GSeq, int startPos, int endPos, uint64_t &MinimizerValue);
uint64_t getSizeofFile(const char *filename);


//void	convertInt64ToString(const uint64_t *GSeqInt, char *GSeq, int char_len);
//int	CheckBadMinSubstring(const char *GSeq, int i);
//int	getMinimizerValue(const char *GSeq, int MinSubPos);
//int	findMinimumPSubstring(const char *GSeq, int startPos, int endPos);

#endif

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <thread>
#include <cstdlib>
#include <memory>

#include <sys/stat.h>
#include <sys/types.h>

#if defined(_WIN32) || defined(_WIN64)
/* We are on Windows */
#include <direct.h>
#include <windows.h>
#include <stdio.h>
#include <tchar.h>

#else

/* We are on Non-Windows */
#include <unistd.h>
#define _rmdir rmdir
#define _mkdir(name) mkdir(name,S_IRWXU)
#endif

#include "mylib.h"

const int BUFFERINCREMENTSIZE = 2000000;
const int MAXLINELENGTH = 256;
const int MAXPARTITION = 256;
const uint64_t MINFILESIZEFORFILTER = 500000000;	// ~500mb
const uint64_t BIGFILESIZE = 10000000000;

int G_MMRLen = 10;	//Global version of mmrLen: minimizer length, might need to change the value


TopKmerCounting::TopKmerCounting(char *filename, const int givenKmerSize, const int givenTopCount)
    :kmersize(givenKmerSize), topcount(givenTopCount), // initializer list for const variable members
    maxLineLenInFile(MAXLINELENGTH), maxPartitionNumber(MAXPARTITION),
    minimizerHistogramDiv(new float[(1 << (G_MMRLen * 2)) + 1]), sortedMinimizersDiv(new float[(1 << (G_MMRLen * 2)) + 1]),
    minimizerHistogramFac(new uint32_t[(1 << (G_MMRLen * 2)) + 1]), sortedMinimizersFac(new uint32_t[(1 << (G_MMRLen * 2)) + 1]),
    minimizerHistogramSum(new uint32_t[(1 << (G_MMRLen * 2)) + 1]), sortedMinimizersSum(new uint32_t[(1 << (G_MMRLen * 2)) + 1])
{
    if (givenKmerSize>90 || givenKmerSize<3)
    {
        std::cerr << "Warning: kmersize value is only allowed in range 3-90" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (givenTopCount<1)
    {
        std::cerr << "Warning: topcount value must be bigger than 0" << std::endl;
        exit(EXIT_FAILURE);
    }

    fastqFilename = filename;
    if (givenKmerSize<G_MMRLen && givenKmerSize != 1) // In case kmersize is less than defaul mmrLen
    {
        G_MMRLen = givenKmerSize - 1;
    }
    if (givenKmerSize == 1)
    {
        G_MMRLen = 1;
    }

    mmrLen = G_MMRLen;

    // 	the parameters below are not really good, with enough time one can get proper
    //	formula for bigfiles when high topcount is asked
    if (getSizeofFile(fastqFilename.c_str()) < MINFILESIZEFORFILTER)
    {
        maxDepthSearch = (1 << (mmrLen * 2)) - 1;			// To be sure toplist is correct
        isBigFileEnabled = 0;
        histogramReadRate = 10;
    }
    else {
        if (topcount<25)
        {
            maxDepthSearch = 50;
        }
        else {
            maxDepthSearch = topcount * 2;
        }
        isBigFileEnabled = 1;
        histogramReadRate = 10;
    }

    if (topcount>1000 && isBigFileEnabled)
    {
        isDiskMethodEnabled = 1;
        maxDepthSearch = topcount * 2;
        histogramReadRate = 10;
    }
    else {
        isDiskMethodEnabled = 0;
    }

    if (topcount>1000 && getSizeofFile(fastqFilename.c_str())>BIGFILESIZE)
    {
        isDiskMethodEnabled = 1;
        maxDepthSearch = topcount * 4;
        histogramReadRate = 5;
    }
    else if (getSizeofFile(fastqFilename.c_str())>BIGFILESIZE)
    {
        isDiskMethodEnabled = 0;
        maxDepthSearch = topcount * 4;
        histogramReadRate = 5;
    }
    // 	the parameters above are not really good, with enough time one can get proper
    //	formula for bigfiles when high topcount is asked. might calculating deviation help?

    if (maxDepthSearch>(1 << (mmrLen * 2)))
    {
        maxDepthSearch = (1 << (mmrLen * 2)) - 1;        //In case maxDepthSearch is bigger than border
    }

    nThreads = std::thread::hardware_concurrency();
    threadFlag = new int[maxPartitionNumber];
    for (int i = 0; i<maxPartitionNumber; i++)
    {
        threadFlag[i] = 0;					//Initially all partitions are not processed so they are initialized 0
    }

    myTopCountTable = new std::multimap<int, std::string>[nThreads];			//all threads need its own sorted list before merging

                                                                                //(1<<(mmrLen*2)) is the number of all different minimizers

    for (int i = 0; i <= (1 << (mmrLen * 2)); i++)
    {
        minimizerHistogramDiv[i] = 0.0;	minimizerHistogramFac[i] = 0; minimizerHistogramSum[i] = 0;
    }

    //filtered data must have been allocated by malloc since it realloc will be used
    filteredData = (char**)malloc(maxPartitionNumber * sizeof(char*));
    for (int i = 0; i<maxPartitionNumber; i++)
    {
        filteredData[i] = (char *)malloc(2 * sizeof(char));    //2 is just a dummy number :)
    }

    currentBufferSize = new uint32_t[maxPartitionNumber];
    currentUsedBufferSize = new uint32_t[maxPartitionNumber];
    for (int i = 0; i<maxPartitionNumber; i++)
    {
        currentBufferSize[i] = 0;
    }
    for (int i = 0; i<maxPartitionNumber; i++)
    {
        currentUsedBufferSize[i] = 0;
    }

}

TopKmerCounting::~TopKmerCounting() {
    myTopCountTable[0].clear();				//	others already cleared when combining

    delete[] myTopCountTable;
    delete[] threadFlag;
    free(filteredData);						//	others already cleared during counting
}

void TopKmerCounting::StartCounting() {
    if (isDiskMethodEnabled)
    {
        RunProcessInDISK();
    }
    else {
        RunProcessInRAM();
    }
}

void TopKmerCounting::RunProcessInRAM() {

    HistogramProcess();					//Histogram function for minimizers

    //std::cout << "hist done" << std::endl;
    std::copy_n(minimizerHistogramDiv.get(), (uint32_t)(1 << (mmrLen * 2)), sortedMinimizersDiv.get());
    std::copy_n(minimizerHistogramFac.get(), (uint32_t)(1 << (mmrLen * 2)), sortedMinimizersFac.get());
    std::copy_n(minimizerHistogramSum.get(), (uint32_t)(1 << (mmrLen * 2)), sortedMinimizersSum.get());

    //We need sorted histogram for later filtering
    std::sort(sortedMinimizersDiv.get(), sortedMinimizersDiv.get() + (uint32_t)(1 << (mmrLen * 2)), [](const float x, const float y)->bool {return x > y;});
    std::sort(sortedMinimizersFac.get(), sortedMinimizersFac.get() + (uint32_t)(1 << (mmrLen * 2)), [](const uint32_t x, const uint32_t y)->bool {return x > y;});
    std::sort(sortedMinimizersSum.get(), sortedMinimizersSum.get() + (uint32_t)(1 << (mmrLen * 2)), [](const uint32_t x, const uint32_t y)->bool {return x > y;});

    //All threads will have initialized toplist map with dummy values
    for (int t = 0; t < nThreads; t++) {
        for (int i = 0; i < topcount; i++)	myTopCountTable[t].insert(std::make_pair(0 - i, "dummyString"));
    }

    //std::cout << "partitionProcess" << std::endl;

    partitionProcess();		// by using histograms, file is written into partitions buffer and also filtered

    //all threads will process different filtered partition data and keep always toplist by inserting into their own hashtable
    std::unique_ptr<std::thread[]> pthrds(new std::thread[nThreads]);
    std::thread *partitionThreads = pthrds.get();
    for (int i = 0; i<nThreads; i++) partitionThreads[i] = std::thread([this, i] { this->partition2Table(i); });
    for (int i = 0; i<nThreads; i++) partitionThreads[i].join();

    // after all threads done, merging toplist maps into myTopCountTable[0]
    for (int i = 1; i<nThreads; i++)
    {
        auto minIt = myTopCountTable[0].begin();
        int myTopCountTableMin = minIt->first;

        for (auto it = myTopCountTable[i].begin(); it != myTopCountTable[i].end(); ++it)
        {
            if (it->first > myTopCountTableMin)
            {
                myTopCountTable[0].erase(minIt);
                myTopCountTable[0].insert(std::make_pair(it->first, it->second));
                minIt = myTopCountTable[0].begin();
                myTopCountTableMin = minIt->first;
            }
        }
        myTopCountTable[i].clear();
    }
}

void TopKmerCounting::RunProcessInDISK() {

    HistogramProcess();					//Histogram function for minimizers

    //std::cout << "hist done" << std::endl;

    std::copy_n(minimizerHistogramDiv.get(), (uint32_t)(1 << (mmrLen * 2)), sortedMinimizersDiv.get());
    std::copy_n(minimizerHistogramFac.get(), (uint32_t)(1 << (mmrLen * 2)), sortedMinimizersFac.get());
    std::copy_n(minimizerHistogramSum.get(), (uint32_t)(1 << (mmrLen * 2)), sortedMinimizersSum.get());

    //We need sorted histogram for later filtering
    std::sort(sortedMinimizersDiv.get(), sortedMinimizersDiv.get() + (uint32_t)(1 << (mmrLen * 2)), [](const float x, const float y)->bool {return x > y;});
    std::sort(sortedMinimizersFac.get(), sortedMinimizersFac.get() + (uint32_t)(1 << (mmrLen * 2)), [](const uint32_t x, const uint32_t y)->bool {return x > y;});
    std::sort(sortedMinimizersSum.get(), sortedMinimizersSum.get() + (uint32_t)(1 << (mmrLen * 2)), [](const uint32_t x, const uint32_t y)->bool {return x > y;});

    //All threads will have initialized toplist map with dummy values
    for (int t = 0; t<nThreads; t++)
    {
        for (int i = 0; i < topcount; i++)
        {
            myTopCountTable[t].insert(std::make_pair(0 - i, "dummyString"));
        }
    }

    //std::cout << "partitionProcess" << std::endl;

    partitionProcessDiskMethod();		// by using histograms, file is written into partitions buffer and also filtered

    //all threads will process different filtered partition data and keep always toplist by inserting into their own hashtable
    std::unique_ptr<std::thread[]> pthrds(new std::thread[nThreads]);
    std::thread *partitionThreads = pthrds.get();
    for (int i = 0; i<nThreads; i++) partitionThreads[i] = std::thread([this,i] { this->partition2TableDiskMethod(i); } );
    for (int i = 0; i<nThreads; i++) partitionThreads[i].join();

    std::string tempDir;
    tempDir.append("./temp");
    remove(tempDir.c_str());

    // after all threads done, merging toplist maps into myTopCountTable[0]
    for (int i = 1; i<nThreads; i++)
    {
        auto minIt = myTopCountTable[0].begin();
        int myTopCountTableMin = minIt->first;

        for (auto it = myTopCountTable[i].begin(); it != myTopCountTable[i].end(); ++it)
        {
            if (it->first > myTopCountTableMin)
            {
                myTopCountTable[0].erase(minIt);
                myTopCountTable[0].insert(std::make_pair(it->first, it->second));
                minIt = myTopCountTable[0].begin();
                myTopCountTableMin = minIt->first;
            }
        }
        myTopCountTable[i].clear();
    }
}

void TopKmerCounting::DisplayTopList() {
    for (auto it = myTopCountTable[0].end(); it != myTopCountTable[0].begin(); )
    {
        it--;
        std::cout << it->second << " " << it->first << std::endl;
    }
}


/**
* Function:  getMinimizerValue(const char *GSeq, int MinSubPos)
* It calculates the integer value of a minimizer 				**
* for given char array and position of minimizer				**
* 																**
* 	**************************************************************
* */

/*
int getMinimizerValue(const char *GSeq, int MinSubPos){
int ret_value = 0;
for(int i=0 ; i<G_MMRLen ; i++){
ret_value *=4;
if(GSeq[MinSubPos+i] == 'C')
ret_value += 1;
else if(GSeq[MinSubPos+i] == 'G')
ret_value += 2;
else if(GSeq[MinSubPos+i] == 'T')
ret_value += 3;
}
return ret_value;
}
*/

/**
int CheckBadMinSubstring(const char *GSeq, int i){
if(!strncmp(GSeq+i,"AAA",3) || !strncmp(GSeq+i,"ACA",3))
return 1;
for(int j=1 ; j<G_MMRLen-1 ; j++)
if(!strncmp(GSeq+i+j,"AA",2))
return 1;
return 0;
}
**/

/**
* Function:	CheckBadMinSubstring(uint64_t )
* 		Int64 version of CheckBadMinSubstring(const char *GSeq, int i).
* 		This function checks if given IntString has AAA and ACA prefixes
*  	and also AA anywhere except beginning of the string
*
* */

static int CheckBadMinSubstring(uint64_t nextCand) {
    if ((nextCand >> ((G_MMRLen - 3) * 2)) == 0 || (nextCand >> ((G_MMRLen - 3) * 2)) == 0x04)
    {        //Checking AAA and ACA prefix
        return 1;
    }
    for (int i = 2; i<(G_MMRLen); i++) {													//Checking AA if it anywhere except beginning
        if ((nextCand & 0xf) == 0)
        {
            return 1;
        }
        nextCand >>= 2;
    }
    return 0;
}


/**
* 	Function:  CompareLastPSubstringWithMin(const uint64_t *, int, int ,uint64_t &, uint64_t &)
*
* 	This function check if the last G_MMRLen length substring
* 	is a new minimizer after shifting string one letter
*
* */

static int CompareLastPSubstringWithMin(const uint64_t *GSeq, int endPos, int MinSubPos, uint64_t &minCand, uint64_t &nextCand) {
    int leftSideBitLen = (((endPos - G_MMRLen) & 0x1f) * 2);
    int intIndex = (endPos - G_MMRLen) >> 5;
    if ((64 - leftSideBitLen) >= (G_MMRLen * 2))
    {
        nextCand = (GSeq[intIndex] >> (64 - leftSideBitLen - (G_MMRLen * 2))) & andTable[(G_MMRLen * 2)];
    }
    else {
        nextCand = (GSeq[intIndex] & andTable[(64 - leftSideBitLen)]);
        nextCand = (nextCand << (((G_MMRLen * 2) + leftSideBitLen) - 64)) ^
            (GSeq[intIndex + 1] >> (128 - ((G_MMRLen * 2) + leftSideBitLen)));
    }
    if (minCand > nextCand)
    {
        if (CheckBadMinSubstring(nextCand))
        {
            return 0;
        }
        return 1;
    }
    return 0;
}

/**
int findMinimumPSubstring(const char *GSeq, int startPos, int endPos){
int MinSubPos = startPos;
for(int i=startPos+1 ; i<=endPos-G_MMRLen ; i++){
for(int j=0 ; j<G_MMRLen ; j++){
if(GSeq[MinSubPos+j] < GSeq[i+j])
break;
else
if(GSeq[MinSubPos+j] > GSeq[i+j]){
if(CheckBadMinSubstring(GSeq, i))
break;
MinSubPos=i;
break;
}
}
}
return MinSubPos;
}
**/

/**
* Function:  findMinimumPSubstring(const uint64_t *GSeq, int startPos, int endPos, uint64_t &MinimizerValue)
* This function is the 64bit integer version of findMinimumPSubstring above.
* for given int64 array, by shifting 2 bits i.e. one letter, all substrings are checked one by one
* faster than char version
* it is a little hard to understand the conversion
* indices are so complex because when a minimizer is in two int64, one needs to shift both int64 correctly
* */

static int findMinimumPSubstring(const uint64_t *GSeq, int startPos, int endPos, uint64_t &MinimizerValue) {
    int MinSubPos = startPos;
    uint64_t minCand, nextCand;
    int leftSideBitLen = ((MinSubPos & 0x1f) * 2);
    int intIndex = (MinSubPos) >> 5;
    if ((64 - leftSideBitLen) >= (G_MMRLen * 2))
    {
        minCand = (GSeq[intIndex] >> (64 - leftSideBitLen - (G_MMRLen * 2))) & andTable[(G_MMRLen * 2)];
    }
    else {
        minCand = GSeq[intIndex] & andTable[(64 - leftSideBitLen)];
        minCand = (minCand << (((G_MMRLen * 2) + leftSideBitLen) - 64)) ^
            (GSeq[intIndex + 1] >> (128 - ((G_MMRLen * 2) + leftSideBitLen)));
    }

    for (int i = startPos + 1; i <= endPos - G_MMRLen; i++) {
        leftSideBitLen = ((i & 0x1f) * 2);
        intIndex = i >> 5;
        if ((64 - leftSideBitLen) >= (G_MMRLen * 2))
        {
            nextCand = (GSeq[intIndex] >> (64 - leftSideBitLen - (G_MMRLen * 2))) & andTable[(G_MMRLen * 2)];
        }
        else {
            nextCand = (GSeq[intIndex] & andTable[(64 - leftSideBitLen)]);
            nextCand = (nextCand << (((G_MMRLen * 2) + leftSideBitLen) - 64)) ^
                (GSeq[intIndex + 1] >> (128 - ((G_MMRLen * 2) + leftSideBitLen)));
        }
        if (minCand > nextCand)
        {
            if (CheckBadMinSubstring(nextCand))
            {
                continue;
            }
            MinSubPos = i;
            minCand = nextCand;
        }
    }
    MinimizerValue = minCand;
    return MinSubPos;
}


/**
* Function:	convertStringToInt64(const char *, uint64_t *)
* This function converts given string into int64 arrays and fills the empty array
* A=00	C=01	G=10	T=11
* */

static void convertStringToInt64(const char *GSeq, uint64_t *GSeqInt) {
    int char_len = (int)strlen(GSeq);
    int intArrayEnd = ((char_len + 31) >> 5);
    for (int i = 0; i<intArrayEnd; i++) GSeqInt[i] = 0;
    for (int i = 0; i<char_len; i++)
    {
        GSeqInt[i >> 5] = (GSeqInt[i >> 5] << 2) + char2IntTable[GSeq[i] - 'A'];
    }

    for (int i = 0; i<(32 - (char_len & 0x1f)); i++)
    {
        GSeqInt[char_len >> 5] = (GSeqInt[char_len >> 5] << 2);
    }
}

/**
* Function:	convertInt64ToString(const uint64_t *, char *, int )
* This functions converts given int64 array into char string
* 00=A	01=C	10=G	11=T
* */

void convertInt64ToString(const uint64_t *GSeqInt, char *GSeq, int char_len) {
    for (int i = 0; i<char_len; i++)
    {
        GSeq[i] = Int2charTable[((GSeqInt[i >> 5] >> (62 - ((i & 0x1f) << 1))) & 0x3)] + 'A';
    }
    GSeq[char_len] = '\0';
}

/**
* Function:	copySkmerToBuffer(char *, int , int , uint64_t &)
* This function copy superkmer into buffer in a custom way
* adding delimiter char '_' between superkmers and at the end there is null character
* eg.
* AAAAACTGCTCGTATTATTCG_ACATCATCTATCACTATCTATCTATC\0
* when a new superkmer added, this is what happens
* AAAAACTGCTCGTATTATTCG_ACATCATCTATCACTATCTATCTATC_TATCTATCTACTTATCT\0
* */

void TopKmerCounting::copySkmerToBuffer(char *myline, int startpos, int endpos, uint64_t &MinimizerValue) {
    uint32_t partNumber = ((uint32_t)MinimizerValue) % this->maxPartitionNumber;
    if (this->currentBufferSize[partNumber] > this->currentUsedBufferSize[partNumber] + endpos - startpos + 4)
    {
        strncpy(this->filteredData[partNumber] + this->currentUsedBufferSize[partNumber], myline + startpos, endpos - startpos + 1);
        this->filteredData[partNumber][this->currentUsedBufferSize[partNumber] + endpos - startpos + 1] = '_';
        this->currentUsedBufferSize[partNumber] = this->currentUsedBufferSize[partNumber] + endpos - startpos + 2;
    }
    else {
        this->filteredData[partNumber] = (char*)realloc(this->filteredData[partNumber], (this->currentBufferSize[partNumber] + BUFFERINCREMENTSIZE) * sizeof(char));
        this->currentBufferSize[partNumber] += BUFFERINCREMENTSIZE;
        strncpy(this->filteredData[partNumber] + this->currentUsedBufferSize[partNumber], myline + startpos, endpos - startpos + 1);
        this->filteredData[partNumber][this->currentUsedBufferSize[partNumber] + endpos - startpos + 1] = '_';
        this->currentUsedBufferSize[partNumber] = this->currentUsedBufferSize[partNumber] + endpos - startpos + 2;
    }
    this->filteredData[partNumber][this->currentUsedBufferSize[partNumber]] = '\0';
}

/**
* This function opens the fastqfile and read one line each time, calculate minimizers and check histograms
* and then if it is in the range, superkmer including that minimizer will be written into buffers
* */

void TopKmerCounting::partitionProcess() {

    std::ifstream MyFile;
    MyFile.open(this->fastqFilename, std::ifstream::in);

    if (!MyFile)
    {
        std::cerr << "Error opening " << this->fastqFilename << std::endl;
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<char[]> shrdmyline(new char[this->maxLineLenInFile]);
    char * myline = shrdmyline.get();
    std::unique_ptr<uint64_t[]> shrdmyIntLine(new uint64_t[((this->maxLineLenInFile + 31) / 32) + 1]);
    uint64_t *myIntLine = shrdmyIntLine.get();

    MyFile.ignore(MAXLINELENGTH, '\n');
    MyFile.getline(myline, this->maxLineLenInFile);

    int line_len;

    uint64_t MinimizerValue, nextCandMin;

    do {
        line_len = (int)strlen(myline);
        if (strspn(myline, "ACGT") != line_len)
        {
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.getline(myline, this->maxLineLenInFile);
            continue;
        }
        convertStringToInt64(myline, myIntLine);
        int min_pos = findMinimumPSubstring(myIntLine, 0, this->kmersize, MinimizerValue);
        int next_min_pos;
        int SKmerPosStart = 0;
        int SKmerPosEnd = this->kmersize - 1;
        for (int i = 1; i<line_len - this->kmersize + 1; i++)
        {
            if (i>min_pos)
            {
                if ((this->minimizerHistogramDiv[((uint32_t)MinimizerValue)] > this->sortedMinimizersDiv[this->maxDepthSearch]) ||
                    (this->minimizerHistogramFac[((uint32_t)MinimizerValue)] > this->sortedMinimizersFac[this->maxDepthSearch]) ||
                    (this->minimizerHistogramSum[((uint32_t)MinimizerValue)] > this->sortedMinimizersSum[this->maxDepthSearch]))
                {
                    next_min_pos = findMinimumPSubstring(myIntLine, i, i + this->kmersize, nextCandMin);
                    if (nextCandMin == MinimizerValue)
                    {
                        min_pos = next_min_pos;
                        SKmerPosEnd++;
                        continue;
                    }
                    copySkmerToBuffer(myline, SKmerPosStart, SKmerPosEnd, MinimizerValue);
                }
                SKmerPosStart = i;
                min_pos = findMinimumPSubstring(myIntLine, i, i + this->kmersize, MinimizerValue);

            }
            else if (CompareLastPSubstringWithMin(myIntLine, i + this->kmersize, min_pos, MinimizerValue, nextCandMin))
            {
                if ((this->minimizerHistogramDiv[((uint32_t)MinimizerValue)] > this->sortedMinimizersDiv[this->maxDepthSearch]) ||
                    (this->minimizerHistogramFac[((uint32_t)MinimizerValue)] > this->sortedMinimizersFac[this->maxDepthSearch]) ||
                    (this->minimizerHistogramSum[((uint32_t)MinimizerValue)] > this->sortedMinimizersSum[this->maxDepthSearch]))
                {
                    copySkmerToBuffer(myline, SKmerPosStart, SKmerPosEnd, MinimizerValue);
                }
                SKmerPosStart = i;
                MinimizerValue = nextCandMin;
                min_pos = i + this->kmersize - this->mmrLen;

            }
            SKmerPosEnd++;
        }
        if ((this->minimizerHistogramDiv[((uint32_t)MinimizerValue)] > this->sortedMinimizersDiv[this->maxDepthSearch]) ||
            (this->minimizerHistogramFac[((uint32_t)MinimizerValue)] > this->sortedMinimizersFac[this->maxDepthSearch]) ||
            (this->minimizerHistogramSum[((uint32_t)MinimizerValue)] > this->sortedMinimizersSum[this->maxDepthSearch]))
        {
            copySkmerToBuffer(myline, SKmerPosStart, SKmerPosEnd, MinimizerValue);
        }
        MyFile.ignore(MAXLINELENGTH, '\n');
        MyFile.ignore(MAXLINELENGTH, '\n');
        MyFile.ignore(MAXLINELENGTH, '\n');
        MyFile.getline(myline, this->maxLineLenInFile);

    } while (MyFile.good());
    MyFile.close();
}


/**
* This function opens the fastqfile and read one line each time, calculate minimizers and check histograms
* and then if it is in the range, superkmer including that minimizer will be written into files
* */

void TopKmerCounting::partitionProcessDiskMethod() {
    char buffer[100];
    std::string tempDir;
    tempDir.append("./temp");
    std::string binFileName("/kmer");
    _rmdir(tempDir.c_str());

    if (_mkdir(tempDir.c_str()) != 0)
    {
        std::cerr << "Permission Denied for MKDIR" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ifstream MyFile;
    std::unique_ptr<std::ofstream[]> shrbinfile(new std::ofstream[this->maxPartitionNumber]);
    std::ofstream *BinFile = shrbinfile.get();
    for (int f = 0; f<this->maxPartitionNumber; f++)
    {
        sprintf(buffer, "%s%s%d.txt", tempDir.c_str(), binFileName.c_str(), f);
        BinFile[f].open(buffer, std::ifstream::trunc);
    }

    MyFile.open(this->fastqFilename, std::ifstream::in);

    if (!MyFile)
    {
        std::cerr << "Error opening " << this->fastqFilename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<char[]> shrmyline(new char[this->maxLineLenInFile]);
    char * myline = shrmyline.get();
    std::unique_ptr<uint64_t[]> shrmyIntLine(new uint64_t[((this->maxLineLenInFile + 31) / 32) + 1]);
    uint64_t *myIntLine = shrmyIntLine.get();

    MyFile.ignore(MAXLINELENGTH, '\n');
    MyFile.getline(myline, this->maxLineLenInFile);

    int line_len;

    uint64_t MinimizerValue, nextCandMin;
    do {
        line_len = (int)strlen(myline);
        if (strspn(myline, "ACGT") != line_len)
        {
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.getline(myline, this->maxLineLenInFile);
            continue;
        }
        convertStringToInt64(myline, myIntLine);
        int min_pos = findMinimumPSubstring(myIntLine, 0, this->kmersize, MinimizerValue);
        int SKmerPosStart = 0;
        int SKmerPosEnd = this->kmersize - 1;
        for (int i = 1; i<line_len - this->kmersize + 1; i++)
        {
            if (i>min_pos)
            {
                if ((this->minimizerHistogramDiv[((uint32_t)MinimizerValue)] > this->sortedMinimizersDiv[this->maxDepthSearch]) ||
                    (this->minimizerHistogramFac[((uint32_t)MinimizerValue)] > this->sortedMinimizersFac[this->maxDepthSearch]) ||
                    (this->minimizerHistogramSum[((uint32_t)MinimizerValue)] > this->sortedMinimizersSum[this->maxDepthSearch]))
                {
                    BinFile[((uint32_t)MinimizerValue) % this->maxPartitionNumber].write(myline + SKmerPosStart, SKmerPosEnd - SKmerPosStart + 1);
                    BinFile[((uint32_t)MinimizerValue) % this->maxPartitionNumber].write("\n", 1);
                }
                SKmerPosStart = i;
                min_pos = findMinimumPSubstring(myIntLine, i, i + this->kmersize, MinimizerValue);

            }
            else if (CompareLastPSubstringWithMin(myIntLine, i + this->kmersize, min_pos, MinimizerValue, nextCandMin))
            {
                if ((this->minimizerHistogramDiv[((uint32_t)MinimizerValue)] > this->sortedMinimizersDiv[this->maxDepthSearch]) ||
                    (this->minimizerHistogramFac[((uint32_t)MinimizerValue)] > this->sortedMinimizersFac[this->maxDepthSearch]) ||
                    (this->minimizerHistogramSum[((uint32_t)MinimizerValue)] > this->sortedMinimizersSum[this->maxDepthSearch]))
                {
                    BinFile[((uint32_t)MinimizerValue) % this->maxPartitionNumber].write(myline + SKmerPosStart, SKmerPosEnd - SKmerPosStart + 1);
                    BinFile[((uint32_t)MinimizerValue) % this->maxPartitionNumber].write("\n", 1);
                }
                SKmerPosStart = i;
                MinimizerValue = nextCandMin;
                min_pos = i + this->kmersize - this->mmrLen;

            }
            SKmerPosEnd++;
        }
        if ((this->minimizerHistogramDiv[((uint32_t)MinimizerValue)] > this->sortedMinimizersDiv[this->maxDepthSearch]) ||
            (this->minimizerHistogramFac[((uint32_t)MinimizerValue)] > this->sortedMinimizersFac[this->maxDepthSearch]) ||
            (this->minimizerHistogramSum[((uint32_t)MinimizerValue)] > this->sortedMinimizersSum[this->maxDepthSearch]))
        {
            BinFile[((uint32_t)MinimizerValue) % this->maxPartitionNumber].write(myline + SKmerPosStart, SKmerPosEnd - SKmerPosStart + 1);
            BinFile[((uint32_t)MinimizerValue) % this->maxPartitionNumber].write("\n", 1);
        }
        MyFile.ignore(MAXLINELENGTH, '\n');
        MyFile.ignore(MAXLINELENGTH, '\n');
        MyFile.ignore(MAXLINELENGTH, '\n');
        MyFile.getline(myline, this->maxLineLenInFile);
    } while (MyFile.good());

    for (int f = 0; f<this->maxPartitionNumber; f++)
    {
        sprintf(buffer, "%s%s%d.txt", tempDir.c_str(), binFileName.c_str(), f);
        BinFile[f].close();
    }
    MyFile.close();
}



/**
* Function:	HashTableProcess(int , int )
* This function run by different threads and each thread process different filtered partition buffer
* gets the superkmer and count each kmer in it and insert it into hashtable or increase its counter
* after that thread will update its own sorted map by checking all elements in hashtable
* */

void TopKmerCounting::HashTableProcess(int partNo, int threadNo) {

    if (this->currentUsedBufferSize[partNo] == 0)
    {
        return;
    }
    std::unordered_map<std::string, int> kmerHashTable;
    kmerHashTable.reserve(this->currentUsedBufferSize[partNo] / this->kmersize);
    std::unique_ptr<char[]> shrkmerRead(new char[this->kmersize + 1]);
    char * kmerRead = shrkmerRead.get();
    char *mySuperkmer, *saveptr;

    mySuperkmer = strtok_r(this->filteredData[partNo], "_", &saveptr);		//strtok_r() is thread safe version of strtok() which is not thread safe
    while (mySuperkmer != NULL)
    {
        for (int j = 0; j<(strlen(mySuperkmer) + 1 - (this->kmersize)) && (this->kmersize) <= strlen(mySuperkmer); j++)
        {
            strncpy(kmerRead, mySuperkmer + j, this->kmersize);
            kmerRead[this->kmersize] = '\0';
            auto it = kmerHashTable.find(kmerRead);
            if (it != kmerHashTable.end()) {
                it->second = it->second + 1;
            }
            else {
                kmerHashTable.insert(std::make_pair(kmerRead, 1));
            }
        }
        mySuperkmer = strtok_r(NULL, "_", &saveptr);							//strtok_r() is thread safe version of strtok() which is not thread safe
    }

    auto minIt = this->myTopCountTable[threadNo].begin();
    int myTopCountTableMin = minIt->first;

    for (auto it = kmerHashTable.begin(); it != kmerHashTable.end(); ++it)
    {
        if (it->second > myTopCountTableMin)
        {
            this->myTopCountTable[threadNo].erase(minIt);
            this->myTopCountTable[threadNo].insert(std::make_pair(it->second, it->first));
            minIt = this->myTopCountTable[threadNo].begin();
            myTopCountTableMin = minIt->first;
        }
    }
    kmerHashTable.clear();

    free(this->filteredData[partNo]);
}


/**
* Function:	HashTableProcessDiskMethod(char *, int )
* Same as HashTableProcess function but reads files instead of buffers
* */
void TopKmerCounting::HashTableProcessDiskMethod(char *Partitionfilename, int threadNo) {
    std::ifstream partitionFile;
    partitionFile.open(Partitionfilename, std::ifstream::in);

    if (!partitionFile.good())
    {
        std::cout << "File open error" << std::endl;
        return;
    }
    if (partitionFile.peek() == std::ifstream::traits_type::eof())
    {
        return;
    }
    std::unordered_map<std::string, int> kmerHashTable;
    std::unique_ptr<char[]> shrmySuperkmer(new char[this->maxLineLenInFile]);
    char * mySuperkmer = shrmySuperkmer.get();
    std::unique_ptr<char[]> shrkmerRead(new char[this->kmersize + 1]);
    char * kmerRead = shrkmerRead.get();

    partitionFile.getline(mySuperkmer, this->maxLineLenInFile, '\n');

    do {
        for (int j = 0; j<(strlen(mySuperkmer) - this->kmersize + 1); j++)
        {
            strncpy(kmerRead, mySuperkmer + j, this->kmersize);
            kmerRead[this->kmersize] = '\0';

            auto it = kmerHashTable.find(kmerRead);
            if (it != kmerHashTable.end())
            {
                it->second = it->second + 1;
            }
            else {
                kmerHashTable.insert(std::make_pair(kmerRead, 1));
            }
        }
        partitionFile.getline(mySuperkmer, this->maxLineLenInFile, '\n');
    } while (partitionFile.good());

    auto minIt = this->myTopCountTable[threadNo].begin();
    int myTopCountTableMin = minIt->first;

    for (auto it = kmerHashTable.begin(); it != kmerHashTable.end(); ++it)
    {
        if (it->second > myTopCountTableMin)
        {
            this->myTopCountTable[threadNo].erase(minIt);
            this->myTopCountTable[threadNo].insert(std::make_pair(it->second, it->first));
            minIt = this->myTopCountTable[threadNo].begin();
            myTopCountTableMin = minIt->first;
        }
    }
    kmerHashTable.clear();
    partitionFile.close();
}



/**
* Function:	partition2Table(int )
* Each thread will run this function and gets the partition buffer which is not processed
* */

void TopKmerCounting::partition2Table(int threadNo) {

    for (int partNo = 0; partNo<this->maxPartitionNumber; partNo++)
    {
        this->partitionMutex.lock();
        if (this->threadFlag[partNo] == 0)
        {
            this->threadFlag[partNo] = 1;
            this->partitionMutex.unlock();
        }
        else {
            this->partitionMutex.unlock();
            continue;
        }
        HashTableProcess(partNo, threadNo);
    }
}


/**
* Function:	partition2TableDiskMethod(TopKmerCounting *,int )
* Each thread will run this function and gets the partition file which is not processed
* */
void TopKmerCounting::partition2TableDiskMethod(int threadNo) {

    char buffer[100];
    std::string tempDir;
    tempDir.append("./temp");
    std::string binFileName("/kmer");

    for (int f = 0; f<this->maxPartitionNumber; f++)
    {
        this->partitionMutex.lock();
        if (this->threadFlag[f] == 0)
        {
            this->threadFlag[f] = 1;
            this->partitionMutex.unlock();
        }
        else {
            this->partitionMutex.unlock();
            continue;
        }
        sprintf(buffer, "%s%s%d.txt", tempDir.c_str(), binFileName.c_str(), f);
        HashTableProcessDiskMethod(buffer, threadNo);
        remove(buffer);
    }
}


/**
* Function:	HistogramProcess()
* This function calculates histograms in two different ways
* 1st: sum of 1/numberofkmers, numberofkmers: which shares the same minimizer in superkmer
* 2nd: sum of numberofkmers,	numberofkmers: which shares the same minimizer in superkmer
*
* some minimizers are unique to one kmer
* and some minimizers are shared by so many kmers
*
* therefore both histogram are usefull to identify top kmers
* */

void TopKmerCounting::HistogramProcess() {

    std::ifstream MyFile;
    MyFile.open(this->fastqFilename, std::ifstream::in);

    if (!MyFile)
    {
        std::cerr << "Error opening " << this->fastqFilename << std::endl;
        exit(EXIT_FAILURE);
    }
    std::unique_ptr<char[]> shrmyline(new char[this->maxLineLenInFile]);
    char * myline = shrmyline.get();
    std::unique_ptr<uint64_t[]> shrmyIntLine(new uint64_t[((this->maxLineLenInFile + 31) / 32) + 1]);
    uint64_t *myIntLine = shrmyIntLine.get();
    MyFile.ignore(MAXLINELENGTH, '\n');
    MyFile.getline(myline, this->maxLineLenInFile);

    int line_len;
    int min_pos;
    uint64_t MinimizerValue, nextCandMin;
    int numberOfKmers;
    do {
        line_len = (int)strlen(myline);
        if (strspn(myline, "ACGT") != line_len)
        {
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.getline(myline, this->maxLineLenInFile);
            continue;
        }
        convertStringToInt64(myline, myIntLine);
        min_pos = findMinimumPSubstring(myIntLine, 0, this->kmersize, MinimizerValue);
        numberOfKmers = 1;
        for (int i = 1; i<line_len - this->kmersize + 1; i++)
        {
            if (i>min_pos)
            {
                this->minimizerHistogramDiv[(uint32_t)MinimizerValue] += (1.0 / ((float)numberOfKmers));
                this->minimizerHistogramFac[(uint32_t)MinimizerValue] += numberOfKmers;
                this->minimizerHistogramSum[(uint32_t)MinimizerValue] ++;
                min_pos = findMinimumPSubstring(myIntLine, i, i + this->kmersize, MinimizerValue);
                numberOfKmers = 1;
            }
            else if (CompareLastPSubstringWithMin(myIntLine, i + this->kmersize, min_pos, MinimizerValue, nextCandMin))
            {
                this->minimizerHistogramDiv[(uint32_t)MinimizerValue] += (1.0 / ((float)numberOfKmers));
                this->minimizerHistogramFac[(uint32_t)MinimizerValue] += numberOfKmers;
                this->minimizerHistogramSum[(uint32_t)MinimizerValue] ++;
                min_pos = i + this->kmersize - this->mmrLen;
                MinimizerValue = nextCandMin;
                numberOfKmers = 1;
            }
            else {
                numberOfKmers++;
            }
        }
        this->minimizerHistogramDiv[(uint32_t)MinimizerValue] += (1.0 / ((float)numberOfKmers));
        this->minimizerHistogramFac[(uint32_t)MinimizerValue] += numberOfKmers;
        this->minimizerHistogramSum[(uint32_t)MinimizerValue] ++;
        for (int j = 0; j<this->histogramReadRate && MyFile.good(); j++)
        {
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.ignore(MAXLINELENGTH, '\n');
            MyFile.getline(myline, this->maxLineLenInFile);
        }
    } while (MyFile.good());
    MyFile.close();
}

uint64_t getSizeofFile(const char *filename) {
    std::ifstream MyFile;
    MyFile.open(filename, std::ifstream::in);
    if (!MyFile)
    {
        std::cerr << "Error opening " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    MyFile.seekg(0, MyFile.end);
    uint64_t length = MyFile.tellg();
    MyFile.close();
    return length;
}

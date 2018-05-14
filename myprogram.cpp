#include <iostream>
#include <cstring>
#include <map>

#include "mylib.h"

int main(int argc, char **argv){
    if(argc < 4){
	std::cerr << "Usage: " << argv[0] << " fastqfilename kmersize topcount" << std::endl;
	return 0;
    }

    TopKmerCounting mykmer(argv[1],atoi(argv[2]),atoi(argv[3]));
    mykmer.StartCounting();
    mykmer.DisplayTopList();
    return 0;

}


#ifndef GETDISTSHPP
#define GETDISTSHPP
#include <stdio.h>
#include <stdlib.h>

#define MAXDEPTH 5
#define	MAXFILENAMELENGTH 200

class gdParams {
public:
	char testName[MAXFILENAMELENGTH], testSpec[MAXFILENAMELENGTH], sourceVCFTemplate[MAXFILENAMELENGTH],
		sourceVCFName[MAXFILENAMELENGTH], plinkName[MAXFILENAMELENGTH], plinkArgs[MAXFILENAMELENGTH], chr[10];
	int startPos, endPos,separation,nSub,addChr,keepVCF,keepPlinkFiles,chunkSize,allowSwap;
	float minMAF;
	int readParms(int argc, char* argv[]);
	int getNextArg(char* nextArg, int argc, char* argv[], FILE* fp[MAXDEPTH], int* depth, int* argNum);
};
#endif


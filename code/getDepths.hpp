#ifndef GETDEPTHSHPP
#define GETDEPTHSHPP
#include <stdio.h>
#include <stdlib.h>

#define MAXDEPTH 5
#define	MAXFILENAMELENGTH 200

class depthParams {
public:
	char testName[MAXFILENAMELENGTH], inputFile[MAXFILENAMELENGTH], sourceVCFTemplate[MAXFILENAMELENGTH];
	int depthPos, addChr, keepVCF;
	int readParms(int argc, char* argv[]);
	int getNextArg(char* nextArg, int argc, char* argv[], FILE* fp[MAXDEPTH], int* depth, int* argNum);
};
#endif


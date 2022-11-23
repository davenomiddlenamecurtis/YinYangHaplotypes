#ifndef FINDYINGYANGHPP
#define FINDYINGYANGHPP
#include <stdio.h>
#include <stdlib.h>

#define MAXDEPTH 5
#define	MAXFILENAMELENGTH 200

class fyyParams {
public:
	char testName[MAXFILENAMELENGTH], distFileName[MAXFILENAMELENGTH];
	int minRun,separation;
	float maxDist;
	char chr[100];
	int readParms(int argc, char* argv[]);
	int getNextArg(char* nextArg, int argc, char* argv[], FILE* fp[MAXDEPTH], int* depth, int* argNum);
};
#endif


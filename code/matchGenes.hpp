#ifndef MATCHGENESHPP
#define MATCHGENESHPP
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>

#define MAXDEPTH 5
#define	MAXFILENAMELENGTH 200

class mgParams {
public:
	char testName[MAXFILENAMELENGTH], inputFile[MAXFILENAMELENGTH], refGenesFile[MAXFILENAMELENGTH], 
		geneRegionsFile[MAXFILENAMELENGTH];
	int readParms(int argc, char* argv[]);
	int getNextArg(char* nextArg, int argc, char* argv[], FILE* fp[MAXDEPTH], int* depth, int* argNum);
};

class mgGene {
public:
	char name[100], chr[3];
	int txStart, txEnd;
	mgGene(const char* n, const char* c, int s, int e);
	mgGene() { ; }
	int input(FILE* fr);
};

extern std::map<std::string, mgGene> geneMap;
#endif


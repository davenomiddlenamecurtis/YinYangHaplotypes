#include "matchGenes.hpp"
#include <string.h>
#include <list>

// g++ -o matchGenes matchGenes.cpp matchGenesUtils.cpp dcerror.cpp -lm

mgParams mgp;
#define MAXLINELENGTH 5000

char line[MAXLINELENGTH];

mgGene::mgGene(const char* n, const char* c, int s, int e)
{
	strcpy(name, n);
	strcpy(chr, c);
	txStart = s;
	txEnd = e;
}

int mgGene::input(FILE* fr)
{
	if (fgets(line, MAXLINELENGTH - 1, fr) && sscanf(line, "%s %d %d %s", chr, &txStart, &txEnd, name) == 4)
		return 1;
	else
		return 0;
}

int readGenes(FILE* fg, std::map<std::string, mgGene> &geneMap)
{
	char chr[100], name[100],*ptr;
	int txStart, txEnd;
	fgets(line, MAXLINELENGTH - 1, fg);
	geneMap.clear();
	while (fgets(line, MAXLINELENGTH - 1, fg))
	{
		sscanf(line, "%*s %*s %s %*s %d %d %*s %*s %*s %*s %*s %*s %s", chr, &txStart, &txEnd, name);
		if (strncmp(chr, "chr", 3))
			continue;
		ptr = chr + 3;
		if ((atoi(ptr) < 1 || atoi(ptr) > 22) && *ptr != 'X')
			continue;
		if (strlen(chr) > 5)
			continue;
		std::map<std::string, mgGene>::iterator queryIter = geneMap.find(name);
		if (queryIter == geneMap.end())
		{
			mgGene g(name, ptr, txStart, txEnd); 
			geneMap[name] = g;
		}
		else
		{
			if (queryIter->second.txStart > txStart)
				queryIter->second.txStart = txStart;
			if (queryIter->second.txEnd > txEnd)
				queryIter->second.txEnd = txEnd;
		}
	}
	return 1;
}

int main(int argc, char* argv[])
{
	float dist;
	char fn[MAXFILENAMELENGTH],chr[3];
	FILE *fi, *fg, *fo,*fr;
	int c,s,e;
	std::map<std::string, mgGene> gm;
	std::map<std::string, mgGene> gma[23];
	std::list<std::string> geneList;
	mgGene g;
	if (!mgp.readParms(argc, argv))
		exit(1);
	if (mgp.refGenesFile[0]) // read genes from refseq file
	{
		fg = fopen(mgp.refGenesFile, "r");
		readGenes(fg, gm);
		fclose(fg);
		fr = fopen(mgp.geneRegionsFile, "w");
		for (std::map<std::string, mgGene>::iterator it = gm.begin(); it != gm.end(); ++it)
		{
			mgGene* g = &it->second;
			fprintf(fr, "%s\t%d\t%d\t%s\n",g->chr,g->txStart,g->txEnd,g->name);
		}
		fclose(fr);
	}

	fr = fopen(mgp.geneRegionsFile, "r");
	while (g.input(fr))
	{
		if (g.chr[0] == 'X')
			c = 22;
		else
			c = atoi(g.chr) - 1;
		gma[c][g.name] = g;
	}
	fclose(fr);
	sprintf(fn, "mg.%s.txt", mgp.testName);
	fo = fopen(fn, "w");
	fi = fopen(mgp.inputFile, "r");
	while (fgets(line, MAXLINELENGTH - 1, fi) && sscanf(line, "%s %d %d", chr, &s,&e)==3)
	{
		geneList.clear();
		if (chr[0] == 'X')
			c = 22;
		else
			c = atoi(chr) - 1;
		for (std::map<std::string, mgGene>::iterator it = gma[c].begin(); it != gma[c].end(); ++it)
		{
			if (it->second.txEnd<s || it->second.txStart>e)
				continue;
			else
				geneList.push_back(it->second.name);
		}
		fprintf(fo, "%s\t%d\t%d\t%d\t", chr, s, e, geneList.size());
		for (std::list<std::string>::iterator gl = geneList.begin(); gl != geneList.end(); ++gl)
			fprintf(fo,"%s ",(*gl).c_str());
		fprintf(fo, "\n");
	}
	fclose(fi);
	fclose(fo);
	return 0;
}

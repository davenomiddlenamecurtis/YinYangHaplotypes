#include "findYinYang.hpp"
#include "dcerror.hpp"
#include <ctype.h>
#include <string.h>
// g++ -o findYinYang findYinYang.cpp findYinYangUtils.cpp dcerror.cpp -lm

// wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_assembly_regions.txt
// file with list of known alternate reference haplotypes and patches

const char *badRegionsFn = "GCA_000001405.29_GRCh38.p14_assembly_regions.txt";
class region {
public:
	int st, en;
	char chr[100],label[100];
};
#define MAXBADREGIONS 1000

region badRegions[MAXBADREGIONS];

fyyParams fyyp;

#define MAXRUN 1000
#define MAXLINELENGTH 1000

int main(int argc, char* argv[])
{
	char line[MAXLINELENGTH],fn[MAXFILENAMELENGTH],*ptr;
	int pos[MAXRUN], length[MAXRUN],runLength,v,vv,nr,r,inRegion;
	float dist[MAXRUN],d;
	FILE *fd, *fo, * fs, * fbs, *fr;
	if (!fyyp.readParms(argc, argv))
		exit(1);
	fr = fopen(badRegionsFn, "r");
	nr = 0;
	while (fgets(line, MAXLINELENGTH, fr))
	{
		if (line[0] == '#')
			continue;
		sscanf(line, "%*s %s %d %d %s",badRegions[nr].chr, &badRegions[nr].st, &badRegions[nr].en, badRegions[nr].label);
		if (!strcmp(badRegions[nr].chr,fyyp.chr))
			++nr;
	}
	if (fyyp.distFileName[0] == '\0' || (fd = fopen(fyyp.distFileName, "r")) == 0)
	{
		dcerror(1, "Need to specify --dist-file-name correctly\n");
		exit(1);
	}
	sprintf(fn, "%s.results.txt", fyyp.testName);
	fo = fopen(fn, "w");
	sprintf(fn, "%s.summ.txt", fyyp.testName);
	fs = fopen(fn, "w");
	sprintf(fn, "%s.badsumm.txt", fyyp.testName);
	fbs = fopen(fn, "w");
	runLength = 0;
	while (fgets(line, MAXLINELENGTH - 1, fd))
	{
		pos[runLength] = atoi(line);
		if (runLength == 0)
			length[runLength] = 0;
		else
			length[runLength] = pos[runLength] - pos[0];
		ptr = line;
		while (*ptr && isspace(*ptr))
			++ptr;
		while (*ptr && !isspace(*ptr))
			++ptr;
		while (*ptr && isspace(*ptr))
			++ptr;
		for (v = 0; v < fyyp.separation; ++v)
		{
			d = atof(ptr);
			if (d<=fyyp.maxDist)
			{
				dist[runLength] = d;
				++runLength;
				break;
			}
			while (*ptr && !isspace(*ptr))
				++ptr;
			while (*ptr && isspace(*ptr))
				++ptr;
		}
		if (v == fyyp.separation) // run over
		{
			if (runLength >= fyyp.minRun)
			{
				inRegion = 0;
				for (r = 0; r < nr; ++r)
				{
					if (pos[0] <= badRegions[r].en && pos[runLength] >= badRegions[r].st)
					{
						inRegion = 1; break;
					}
				}
				dist[runLength] = runLength;
				if (inRegion)
				{
					fprintf(fbs, "%s\t%d\t%d\t%d\t%d\t%s\n", fyyp.chr, pos[0], pos[runLength], (int)dist[runLength], length[runLength], badRegions[r].label);
				}
				else 
				{
					fprintf(fs, "%s\t%d\t%d\t%d\t%d\n", fyyp.chr, pos[0], pos[runLength], (int)dist[runLength], length[runLength]);
					++runLength;
					for (v = 0; v < runLength; ++v)
						fprintf(fo, "%s\t%d\t%6.4f\t%d\n", fyyp.chr, pos[v], dist[v], length[v]);
				}
			}
			runLength = 0;
		}
		else for (vv = 0; vv < v; ++vv)
			fgets(line, MAXLINELENGTH - 1, fd);
	}
	fclose(fo);
	fclose(fs);
	fclose(fbs);
	fclose(fd);
}

#include "getDepths.hpp"
#include <ctype.h>
#include <string.h>
#include <math.h>

#ifdef TRACE
#define trace printf("Got to line %d OK\n",__LINE__)
#else
#define trace
#endif

// g++ -o getDepths getDepths.cpp getDepthsUtils.cpp dcerror.cpp -lm

// read in a list of coordinates and output average read depth, allele balance, etc

// also not that this will not exclude children in 1000G

depthParams dp;

#define MAXLINELENGTH 200000
// this is enough for 1000G
// for bigger datasets would need to read calls sequentially from file rather than use
// line-based input


int main(int argc, char* argv[])
{
	char *genoLine, fn[MAXFILENAMELENGTH], line[MAXFILENAMELENGTH], commLine[MAXFILENAMELENGTH], 
		chr[3],tempVcf[MAXFILENAMELENGTH], *testAlt,
		*ptr;
	int nSub[3], pos, testPos, i, dd, d[2], f,g,gg;
	float td[3][2], td2[3][2], ab, tab[3], tab2[3];
	FILE *fi, *fv, *fo;
	if (!dp.readParms(argc, argv))
		exit(1);
	genoLine = (char*)malloc(MAXLINELENGTH);
	testAlt = (char*)malloc(MAXLINELENGTH); // there may be a long insertion
	fi = fopen(dp.inputFile, "r");
	sprintf(fn, "results.%s.txt", dp.testName);
	fo = fopen(fn, "w");
	while (fgets(line, MAXFILENAMELENGTH - 1, fi) && sscanf(line, "%s %d", chr, &pos)==2)
	{
		strcpy(fn, dp.sourceVCFTemplate);
		*strstr(fn, "CHR") = '\0';
		sprintf(strchr(fn, '\0'), "%s%s", chr, strstr(dp.sourceVCFTemplate, "CHR") + 3);
		sprintf(tempVcf, "temp.%s.%d.vcf", chr, pos);
		sprintf(commLine, "tabix %s %s%s:%d-%d > %s", fn, dp.addChr ? "chr" : "", chr, pos, pos, tempVcf);
		system(commLine);
		trace;
		fv = fopen(tempVcf, "r");
		trace;
		while (fgets(genoLine, MAXLINELENGTH - 1, fv))
		{
			trace;
			sscanf(genoLine, "%*s %d %*s %*s %s",&testPos,testAlt);
			trace;
			if (testPos != pos || strchr(testAlt, ','))
				continue;
			trace;
			ptr = genoLine;
			for(gg = 0; gg < 3; ++gg)
			{
				nSub[gg] = tab[gg] = tab2[gg] = 0;
				for (dd = 0; dd < 2; ++dd)
					td[gg][dd] = td2[gg][dd] = 0;
			}
			trace;
			for (i=0;i<9;++i)
				{
				while (!isspace(*ptr)) ++ptr;
				while (isspace(*ptr)) ++ptr;
				}
			trace;
			while (*ptr)
			{
				trace;
				if (*ptr != '0' && *ptr != '1') // unknown genotype
				{
					while (!isspace(*ptr)) ++ptr;
					trace;
				}
				else 
				{
					trace;
					if (ptr[1]==':') // male on X chromosome
						g = (ptr[0] - '0')*2; // count with homozygotes
					else
						g = ptr[0] - '0' + ptr[2] - '0';
					++nSub[g];
					for (f = 0; f < dp.depthPos; ++f)
					{
						while (*ptr != ':') ++ptr;
						++ptr;
					}
					sscanf(ptr, "%d,%d", &d[0], &d[1]);
					for (dd = 0; dd < 2; ++dd)
					{
						td[g][dd] += d[dd];
						td2[g][dd] += d[dd] * d[dd];
					}
					ab = float(d[1]) / (d[0] + d[1]);
					tab[g] += ab;
					tab2[g] += ab * ab;
					while (*ptr && !isspace(*ptr)) ++ptr;
				}
				trace;
				while (*ptr && isspace(*ptr)) ++ptr;
				trace;
			}
			break; //  because we found a valid variant and do not need to keep looking
		}
		fclose(fv);
		if (!dp.keepVCF)
			remove(tempVcf);
		for (gg = 0; gg < 3; ++gg)
		{
			for (dd = 0; dd < 2; ++dd)
			{
				td[gg][dd] = td[gg][dd] / nSub[gg]; // mean
				td2[gg][dd] = td2[gg][dd] / nSub[gg] - td[gg][dd] * td[gg][dd];
			}
			tab[gg] = tab[gg] / nSub[gg];
			tab2[gg] = tab2[gg] / nSub[gg] - tab[gg] * tab[gg];
		}
		fprintf(fo, "%s\t%d\t%d\t%d\t%d\t", chr, pos, nSub[0], nSub[1], nSub[2]);
		for (gg = 0; gg < 3; ++gg) 
		{
			if (nSub[gg] == 0)
				fprintf(fo, "0.0\t0.0\t0.00\t0.00\t0.00\t0.00\t");
			else
				fprintf(fo, "%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t",
					td[gg][0], td[gg][1], sqrt(td2[gg][0]), sqrt(td2[gg][1]), tab[gg], sqrt(tab2[gg]));
		}
		fprintf(fo, "\n");
	}
	fclose(fi);
	fclose(fo);
}

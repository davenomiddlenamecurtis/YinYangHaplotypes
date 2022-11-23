#include "findYinYang.hpp"
#include <ctype.h>
// g++ - o findYinYang findYinYang.cpp findYinYangUtils.cpp dcerror.cpp - lm

fyyParams fyyp;

int main(int argc, char* argv[])
{
	char **genoLines,**linePtr,*keptPtr,commLine[1000],fn[MAXFILENAMELENGTH], * ptr1, * ptr2;
	int i, chunkStart, chunkEnd, ch, linesFilled, maxLineLength, thisPos, r, used;
	float dist;
	FILE *fc, * fo;
	if (!fyyp.readParms(argc, argv))
		exit(1);
	genoLines = (char**)calloc(fyyp.separation + 1, sizeof(char*));
	linePtr = (char**)calloc(fyyp.separation + 1, sizeof(char*));
	linesFilled = 0;
	sprintf(fn, "results.%s.txt", fyyp.testSpec);
	fo = fopen(fn, "w");
	for (chunkStart = fyyp.startPos; chunkStart < fyyp.endPos; chunkStart += fyyp.chunkSize)
	{
		chunkEnd = chunkStart + fyyp.chunkSize - 1;
		sprintf(commLine, "tabix -h %s %s%s:%d-%d > %s.%s.%d.%d.vcf",
			fyyp.sourceVCFName, fyyp.addChr ? "chr" : "", fyyp.chr, chunkStart, chunkEnd, 
			fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
// I could add dummy rs names if needed
		puts(commLine);
		system(commLine);
//		sprintf(commLine, "%s %s --vcf %s.%s.%d.%d.vcf --maf %f --recode A-transpose --out %s.%s.%d.%d",
// remove does not work correctly so do this in two stages
		sprintf(commLine, "%s %s --vcf %s.%s.%d.%d.vcf --maf %f --make-bed --out temp",
			fyyp.plinkName, fyyp.plinkArgs,
			fyyp.testName, fyyp.chr, chunkStart, chunkEnd,
			fyyp.minMAF);
		puts(commLine);
		system(commLine);
		sprintf(commLine, "%s --bfile temp --recode A-transpose --out %s.%s.%d.%d",
			fyyp.plinkName, 
			fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
		puts(commLine);
		system(commLine);
		sprintf(fn, "%s.%s.%d.%d.nosex", fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
		fc = fopen(fn, "r");
		for (fyyp.nSub = 0; fgets(commLine, 200, fc); ++fyyp.nSub);
		fclose(fc);
		sprintf(fn, "%s.%s.%d.%d.traw", fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
		fc = fopen(fn, "r");
		// skip first line, which may be long
		do {
			ch = fgetc(fc);
		} while (ch != '\n' && ch != EOF);
		ch = fgetc(fc);
		if (ch == EOF || ch == '\n')
			break; // no more variants in VCF
		// this may fail if first chunk is too near start of chr to catch any variants
		ungetc(ch, fc);
		while (1)
		{
			if (!linesFilled)
			{
				maxLineLength = fyyp.nSub * 3 + 100; // multiply by 3 because genotype can be NA
				for (i = 0; i < fyyp.separation + 1; ++i)
					genoLines[i] = (char*)calloc(maxLineLength, sizeof(char));
				for (i = 0; i < fyyp.separation; ++i)
					fgets(genoLines[i + 1], maxLineLength, fc);
				for (i = 0; i < fyyp.separation + 1; ++i)
					linePtr[i] = genoLines[i];
				linesFilled = 1;
			}
			keptPtr = linePtr[0];
			if (!fgets(linePtr[0], maxLineLength, fc))
				// this chunk used up
			{
				fclose(fc);
				if (!fyyp.keepVCF)
				{
					sprintf(fn, "%s.%s.%d.%d.vcf", fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
					remove(fn);
				}
				if (!fyyp.keepPlinkFiles)
				{
					sprintf(fn, "%s.%s.%d.%d.nosex", fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
					remove(fn);
					sprintf(fn, "%s.%s.%d.%d.traw", fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
					remove(fn);
					sprintf(fn, "%s.%s.%d.%d.log", fyyp.testName, fyyp.chr, chunkStart, chunkEnd);
					remove(fn);
				}
				break;
			}
			for (i = 0; i < fyyp.separation; ++i)
				linePtr[i] = linePtr[i + 1];
			linePtr[fyyp.separation] = keptPtr;
			ptr1 = linePtr[0];
			sscanf(ptr1, "%*s %*s %*s %d", &thisPos);
			fprintf(fo, "%d\t", thisPos);
			for (int v = 0; v < fyyp.separation; ++v)
			{
				ptr1 = linePtr[0];
				ptr2 = linePtr[1+v];
				used = 0;
				dist = 0;
				for (i = 0; i < 6; ++i)
				{
					while (!isspace(*ptr1)) ++ptr1;
					while (isspace(*ptr1)) ++ptr1;
				}
				for (i = 0; i < 6; ++i)
				{
					while (!isspace(*ptr2)) ++ptr2;
					while (isspace(*ptr2)) ++ptr2;
				}
				for (i = 0; i < fyyp.nSub; ++i)
				{
					char g1, g2;
					g1 = *ptr1;
					g2 = *ptr2;
					if (g1 >= '0' && g1 <= '2' && g2 >= '0' && g2 <= '2')
					{
						dist += (g1 > g2) ? g1 - g2 : g2 - g1;
						++used;
					}
					while (!isspace(*ptr1)) ++ptr1; ++ptr1; // can be NA
					while (!isspace(*ptr2)) ++ptr2; ++ptr2;
				}
				fprintf(fo,"%5.3f\t", dist / used);
			}
			fprintf(fo, "\n");
		}
	}
}

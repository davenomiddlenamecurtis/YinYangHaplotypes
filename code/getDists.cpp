#include "getDists.hpp"
#include <ctype.h>
// g++ -o getDists getDists.cpp getDistsUtils.cpp dcerror.cpp -lm

gdParams gdp;

int main(int argc, char* argv[])
{
	char **genoLines,**linePtr,*keptPtr,commLine[1000],fn[MAXFILENAMELENGTH], * ptr1, * ptr2;
	int i, chunkStart, chunkEnd, ch, linesFilled, maxLineLength, thisPos, r, used;
	float dist,hetDist,sameDist,swappedDist;
	FILE *fc, * fo;
	if (!gdp.readParms(argc, argv))
		exit(1);
	genoLines = (char**)calloc(gdp.separation + 1, sizeof(char*));
	linePtr = (char**)calloc(gdp.separation + 1, sizeof(char*));
	linesFilled = 0;
	sprintf(fn, "results.%s.txt", gdp.testSpec);
	fo = fopen(fn, "w");
	for (chunkStart = gdp.startPos; chunkStart < gdp.endPos; chunkStart += gdp.chunkSize)
	{
		chunkEnd = chunkStart + gdp.chunkSize - 1;
		sprintf(commLine, "tabix -h %s %s%s:%d-%d > %s.%s.%d.%d.vcf",
			gdp.sourceVCFName, gdp.addChr ? "chr" : "", gdp.chr, chunkStart, chunkEnd, 
			gdp.testName, gdp.chr, chunkStart, chunkEnd);
// I could add dummy rs names if needed
		fprintf(stderr, "%s\n", commLine); fflush(stderr);
		system(commLine);
//		sprintf(commLine, "%s %s --vcf %s.%s.%d.%d.vcf --maf %f --recode A-transpose --out %s.%s.%d.%d",
// remove does not work correctly so do this in two stages
		sprintf(commLine, "%s %s --vcf %s.%s.%d.%d.vcf --maf %f --make-bed --out temp.%s.%d.%d",
			gdp.plinkName, gdp.plinkArgs,
			gdp.testName, gdp.chr, chunkStart, chunkEnd,
			gdp.minMAF, gdp.chr, chunkStart, chunkEnd);
		fprintf(stderr, "%s\n", commLine); fflush(stderr);
		system(commLine);
		sprintf(fn, "temp.%s.%d.%d.bed", gdp.chr, chunkStart, chunkEnd);
		fc = fopen(fn,"rb");
		if (fc == NULL)
		{
			fprintf(stderr, "%s was not written\n", fn); fflush(stderr);
			if (!gdp.keepVCF)
			{
				sprintf(fn, "%s.%s.%d.%d.vcf", gdp.testName, gdp.chr, chunkStart, chunkEnd);
				remove(fn);
				sprintf(commLine, "rm temp.%s.%d.%d*", gdp.chr, chunkStart, chunkEnd);
				fprintf(stderr, "%s\n", commLine); fflush(stderr);
				system(commLine);
			}
			continue; // vcf file has no variants. would cause fatal error for below plink command
		}
		else
			fclose(fc);
		sprintf(commLine, "%s --bfile temp.%s.%d.%d --recode A-transpose --out %s.%s.%d.%d",
			gdp.plinkName, gdp.chr, chunkStart, chunkEnd,
			gdp.testName, gdp.chr, chunkStart, chunkEnd);
		fprintf(stderr, "%s\n", commLine); fflush(stderr);
		system(commLine);
		sprintf(commLine, "rm temp.%s.%d.%d*", gdp.chr, chunkStart, chunkEnd);
		// just assume unix - the files are quite big
		fprintf(stderr, "%s\n", commLine); fflush(stderr);
		system(commLine);
		sprintf(fn, "%s.%s.%d.%d.nosex", gdp.testName, gdp.chr, chunkStart, chunkEnd);
		fc = fopen(fn, "r");
		for (gdp.nSub = 0; fgets(commLine, 200, fc); ++gdp.nSub);
		fclose(fc);
		sprintf(fn, "%s.%s.%d.%d.traw", gdp.testName, gdp.chr, chunkStart, chunkEnd);
		fc = fopen(fn, "r");
		// skip first line, which may be long
		do {
			ch = fgetc(fc);
		} while (ch != '\n' && ch != EOF);
		ch = fgetc(fc);
		if (ch == EOF || ch == '\n')
			continue; // no variants in VCF
		// if break this may fail if first chunk is too near start of chr to catch any variants
		ungetc(ch, fc);
		while (1)
		{
			if (!linesFilled)
			{
				maxLineLength = gdp.nSub * 3 + 100; // multiply by 3 because genotype can be NA
				for (i = 0; i < gdp.separation + 1; ++i)
					genoLines[i] = (char*)calloc(maxLineLength, sizeof(char));
				for (i = 0; i < gdp.separation; ++i)
					fgets(genoLines[i + 1], maxLineLength, fc);
				for (i = 0; i < gdp.separation + 1; ++i)
					linePtr[i] = genoLines[i];
				linesFilled = 1;
			}
			keptPtr = linePtr[0];
			if (!fgets(linePtr[0], maxLineLength, fc))
				// this chunk used up
			{
				fclose(fc);
				if (!gdp.keepVCF)
				{
					sprintf(fn, "%s.%s.%d.%d.vcf", gdp.testName, gdp.chr, chunkStart, chunkEnd);
					remove(fn);
				}
				if (!gdp.keepPlinkFiles)
				{
					sprintf(fn, "%s.%s.%d.%d.nosex", gdp.testName, gdp.chr, chunkStart, chunkEnd);
					remove(fn);
					sprintf(fn, "%s.%s.%d.%d.traw", gdp.testName, gdp.chr, chunkStart, chunkEnd);
					remove(fn);
					sprintf(fn, "%s.%s.%d.%d.log", gdp.testName, gdp.chr, chunkStart, chunkEnd);
					remove(fn);
				}
				break;
			}
			for (i = 0; i < gdp.separation; ++i)
				linePtr[i] = linePtr[i + 1];
			linePtr[gdp.separation] = keptPtr;
			ptr1 = linePtr[0];
			sscanf(ptr1, "%*s %*s %*s %d", &thisPos);
			fprintf(fo, "%d\t", thisPos);
			for (int v = 0; v < gdp.separation; ++v)
			{
				ptr1 = linePtr[0];
				ptr2 = linePtr[1+v];
				used = 0;
				// dist = 0;
				hetDist = sameDist = swappedDist = 0;
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
				for (i = 0; i < gdp.nSub; ++i)
				{
					char g1, g2;
					g1 = *ptr1;
					g2 = *ptr2;
					if (g1 >= '0' && g1 <= '2' && g2 >= '0' && g2 <= '2')
					{
						// dist += (g1 > g2) ? g1 - g2 : g2 - g1;
						++used;
						if (g1=='1' || g2=='1')
						{ 
							if (g1 != g2)
								++hetDist; // one hom, one het
						}
						else // both hom
						{
							if (g1 == g2)
								++swappedDist;
							else
								++sameDist;
						}
					}
					while (!isspace(*ptr1)) ++ptr1; ++ptr1; // can be NA
					while (!isspace(*ptr2)) ++ptr2; ++ptr2;
				}
				if (gdp.allowSwap)
					dist = hetDist + 2*(swappedDist > sameDist ? sameDist : swappedDist);
				else
					dist = hetDist + 2 * sameDist;
				fprintf(fo,"%7.5f\t", dist / used);
			}
			fprintf(fo, "\n");
		}
	}
}

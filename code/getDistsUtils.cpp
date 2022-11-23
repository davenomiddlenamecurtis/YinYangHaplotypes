
#include "getDists.hpp"
#include "dcerror.hpp"
#include <ctype.h>
#include <string.h>

#define isArgType(a) (a[0]=='-' && a[1]=='-')
// #define FILLARG(str) (strcmp(arg,str) ? 0 : ((getNextArg(arg, argc, argv, fp,&depth, &argNum) && !isArgType(arg)) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))
// do not check if arg type because need to be able to provide arguments for plink
#define FILLARG(str) (strcmp(arg,str) ? 0 : (getNextArg(arg, argc, argv, fp,&depth, &argNum) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))

int gdParams::readParms(int argc, char* argv[])
{
	FILE* fp[MAXDEPTH];
	int argNum,depth;
	char arg[2000];
	depth = -1;
	argNum = 1;
	strcpy(testName, "distances");
	strcpy(sourceVCFTemplate, "");
	strcpy(plinkName, "plink");
	strcpy(plinkArgs, "--remove children.txt");
	strcpy(chr, "");
	startPos=endPos=separation=nSub=addChr=keepVCF=keepPlinkFiles=chunkSize=allowSwap=0;
	minMAF = 0;
	while (getNextArg(arg, argc, argv, fp, &depth, &argNum))
	{
		if (!isArgType(arg))
		{
			dcerror(1, "Expected argument beginning -- but got this: %s\n", arg);
			return 0;
		}
		else if (FILLARG("--arg-file"))
		{
			if (++depth >= MAXDEPTH)
			{
				dcerror(1, "Attempting to recurse too deeply into arg-files with this one: %s\n", arg);
				return 0;
			}
			else
			{
				fp[depth] = fopen(arg, "r");
				if (fp[depth] == NULL)
				{
					dcerror(1, "Could not open arg file: %s\n", arg);
					return 0;
				}
			}
		}
		else if (FILLARG("--test-name"))
		{
			strcpy(testName, arg);
		}
		else if (FILLARG("--vcf-spec"))
		{
			strcpy(sourceVCFTemplate, arg);
		}
		else if (FILLARG("--chr"))
		{
			strcpy(chr, arg);
		}
		else if (FILLARG("--plink"))
		{
			strcpy(plinkName, arg);
		}
		else if (FILLARG("--plink-args"))
		{
			strcpy(plinkArgs, arg);
		}
		else if (FILLARG("--allow-swap"))
		{
			allowSwap = atoi(arg);
		}
		else if (FILLARG("--start-pos"))
		{
			startPos = atoi(arg);
		}
		else if (FILLARG("--end-pos"))
		{
			endPos = atoi(arg);
		}
//		else if (FILLARG("--num-subs"))
//		{
//			nSub = atoi(arg);
//		}
		else if (FILLARG("--separation"))
		{
			separation = atoi(arg);
		}
		else if (FILLARG("--add-chr"))
		{
			addChr = atoi(arg);
		}
		else if (FILLARG("--chunk-size"))
		{
			chunkSize = atoi(arg);
		}
		else if (FILLARG("--keep-vcf"))
		{
			keepVCF = atoi(arg);
		}
		else if (FILLARG("--keep-plink-files"))
		{
			keepPlinkFiles = atoi(arg);
		}
		else if (FILLARG("--min-maf"))
		{
			minMAF = atof(arg);
		}
		else
		{
			dcerror(1, "Unrecognised argument: %s\n", arg);
			return 0;
		}
	}
	strcpy(sourceVCFName, sourceVCFTemplate);
	*strstr(sourceVCFName, "CHR") = '\0';
	sprintf(strchr(sourceVCFName, '\0'), "%s%s", chr, strstr(sourceVCFTemplate, "CHR") + 3);
	sprintf(testSpec, "%s.%s.%d.%d", testName, chr, startPos, endPos);
	return 1;
}

int gdParams::getNextArg(char* nextArg, int argc, char* argv[], FILE* fp[MAXDEPTH], int* depth, int* argNum)
{
	char* ptr;
	int ch;
	*nextArg = '\0';
	while (*depth > -1)
	{
		do {
			ch = fgetc(fp[*depth]);
		} while (ch != EOF && isspace(ch));
		if (ch == '\'')
		{
			ptr = nextArg;
			while ((ch = fgetc(fp[*depth])) != '\'' && ch != EOF)
				*ptr++ = ch;
			*ptr = '\0';
			return 1;
		}
		else if (ch != EOF)
		{
			nextArg[0] = ch;
			ptr = nextArg + 1;
			while ((ch = fgetc(fp[*depth])) != EOF && !isspace(ch))
				*ptr++ = ch;
			*ptr = '\0';
			return 1;
		}
		else
		{
			fclose(fp[*depth]);
			--* depth;
		}
	}
	if (*argNum < argc)
	{
		if (argv[*argNum][0] == '\'' && (ptr = strchr(argv[*argNum] + 1, '\'')) != 0)
		{
			*ptr = '\0';
			strcpy(nextArg, argv[*argNum] + 1);
		}
		else
			strcpy(nextArg, argv[*argNum]);
		++* argNum;
		return 1;
	}
	else
		return 0;
}

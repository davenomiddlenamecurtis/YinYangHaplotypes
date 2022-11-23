
#include "findYinYang.hpp"
#include "dcerror.hpp"
#include <string.h>
#include <ctype.h>
#include <math.h>

#define isArgType(a) (a[0]=='-' && a[1]=='-')
// #define FILLARG(str) (strcmp(arg,str) ? 0 : ((getNextArg(arg, argc, argv, fp,&depth, &argNum) && !isArgType(arg)) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))
// do not check if arg type because need to be able to provide arguments for plink
#define FILLARG(str) (strcmp(arg,str) ? 0 : (getNextArg(arg, argc, argv, fp,&depth, &argNum) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))

int fyyParams::readParms(int argc, char* argv[])
{
	FILE* fp[MAXDEPTH];
	int argNum,depth;
	char arg[2000];
	chr[0] = '\0';
	depth = -1;
	argNum = 1;
	strcpy(testName, "fyy");
	strcpy(distFileName, "");
	minRun = 1000;
	maxDist = 0;
	separation = 0;
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
		else if (FILLARG("--chr"))
		{
			strcpy(chr, arg);
		}
		else if (FILLARG("--test-name"))
		{
			strcpy(testName, arg);
		}
		else if (FILLARG("--dist-file-name"))
		{
			strcpy(distFileName, arg);
		}
		else if (FILLARG("--min-run"))
		{
			minRun = atoi(arg);
		}
		else if (FILLARG("--max-dist"))
		{
			maxDist = atof(arg);
		}
		else if (FILLARG("--separation"))
		{
			separation = atoi(arg);
		}
		else
		{
			dcerror(1, "Unrecognised argument: %s\n", arg);
			return 0;
		}
	}
	if (separation <= 0)
	{
		dcerror(1, "Must supply a value for --separation\n");
		return 0;
	}
	if (chr[0]=='\0')
	{
		dcerror(1, "Must supply a value for --chr\n");
		return 0;
	}
	return 1;
}

int fyyParams::getNextArg(char* nextArg, int argc, char* argv[], FILE* fp[MAXDEPTH], int* depth, int* argNum)
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

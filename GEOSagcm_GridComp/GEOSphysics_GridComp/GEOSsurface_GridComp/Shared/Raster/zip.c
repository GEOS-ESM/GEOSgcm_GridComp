/*  $Id$  */

#include <stdio.h>
#include "zlib.h"

#ifdef STDC
#  include <string.h>
#  include <stdlib.h>
#endif

#define MAXREC 1024

gzFile   file = NULL;

void ZIPBUFF(char *fname, char *Buff, int *nbytes, int *form, int *Write)
{

  char mode[]=" b6";

  mode[0]=Write?'w':'r';

  if(*nbytes) {  // Do I/O
    if(file==NULL) {
      fprintf(stderr, "error: gz file not open %s\n",fname);
      exit(2);
    }

    if(*form) {
      if(Write) {
	(void)gzputs(file,Buff);
	(void)gzputc(file,'\n');
      }
      else
	(void)gzgets(file,Buff,MAXREC);
    }
    else {
      if(Write) 
	(void)gzwrite(file,Buff,*nbytes);
      else
	(void)gzread(file,Buff,*nbytes);
    }

  } else {  // Toggle Open or Close

    if(file==NULL) { 
      if ((file = gzopen(fname,mode)) == NULL) {
	fprintf(stderr, "error: could not open gz file %s\n",fname);
	exit(1);
      }

    } else {
      gzclose(file);
      file=NULL;
    }

  }
}


void zipbuff(char *fname, char *Buff, int *nbytes, int *form, int *Write)
{
  void ZIPBUFF(char *fname, char *Buff, int *nbytes, int *form, int *Write);
}
void zipbuff_(char *fname, char *Buff, int *nbytes, int *form, int *Write)
{
  void ZIPBUFF(char *fname, char *Buff, int *nbytes, int *form, int *Write);
}
void ZIPBUFF_(char *fname, char *Buff, int *nbytes, int *form, int *Write)
{
  void ZIPBUFF(char *fname, char *Buff, int *nbytes, int *form, int *Write);
}

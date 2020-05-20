/*
 * $Header: /u/gcmpack/MITgcm_contrib/PRM/multi_comp_setup/comp_mitgcm/code_basic/setdir.c,v 1.1.1.1 2006/10/10 18:17:28 cnh Exp $
 * $Name:  $
*/

#include "FC_NAMEMANGLE.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

static char popto[4096] = "\0";
                                                                                
void FC_NAMEMANGLE(mysetdir) (char *theDir )
{
   char RUN_DIR[1024];
   char *rundirsetting;
   char *cwdrc;
   char *SROOT;
   char *SPREF;
   char *scycle;
   char *endPtr;
   char SROOT_BUFFER[1024];
   int  rc;
   int  sDirNum;
   int sDirCycle;

   cwdrc = getcwd(popto, 4096);
   sprintf(RUN_DIR,"%s",theDir);
   rc = chdir(RUN_DIR);
     
   if ( rc != 0 ) {
    fprintf(stderr,"\"%s\"\n",theDir);
    //    perror(strerror(rc));
    exit(-1);
   }
     
}
void FC_NAMEMANGLE(popdir)()
{
   char RUN_DIR[1024];
   char *rundirsetting;
   char *cwdrc;
   char *SROOT;
   char *SPREF;
   char *scycle;
   char *endPtr;
   char SROOT_BUFFER[1024];
   int  rc;
   int  sDirNum;
   int sDirCycle;

   if ( popto[0] == '\0' ) {
    return;
   }
   rc = chdir(popto);
}

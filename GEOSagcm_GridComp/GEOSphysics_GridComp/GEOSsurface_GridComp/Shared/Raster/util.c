#include "stdio.h"
#include "stdlib.h"
#include "errno.h"
#include "stdarg.h"
#include "fcntl.h"
#include "string.h"

int Open=0,Form=0,Close=0,Write=1,Read=0;

void WRITERST(int *rst, int *nx, int *ny, char *filename, int *zip)
{
  int fd, i, j, *buf, nbytes;
  int err,  reclen;

  nbytes=*nx*4;
  reclen=nbytes+8;

  if(!(buf = (int *)malloc(reclen))){
    printf("WRTRST: Could not allocate buffer while writing file %s\n",filename);
    exit(1);
  }

  buf[    0]=nbytes;
  buf[*nx+1]=nbytes;

  if(*zip) {
    (void)ZIPBUFF(filename,rst,&Open,&Form,&Write);
  }
  else {
    if((fd=creat(filename,0755))<0) {
      printf("WRTRST: Could not open %s\n",filename);
      exit(1);
    }
  }

  for (j=0;j<*ny;j++) {
    memcpy((char *)(&(buf[1])),(char *)rst,nbytes);
    rst+=*nx;

    if(*zip) {
      (void)ZIPBUFF(filename,buf,&reclen,&Form,&Write);
    }
    else {
      if((err=write(fd,buf,nbytes+8))<0) {
	printf("WRTRST: Write error %d on file %s\n",err,filename);
	exit(1);
      }
    }
  }

  if(*zip) {
    (void)ZIPBUFF(filename,rst,&Close,&Form,&Write);
  }
  else {
    close(fd);
  }

  free(buf);
}


void READRST(int *rst, int *nx, int *ny, char *filename, int *zip)
{
  int j;
  int fd, nbytes, *buf, reclen, intbuf[2];

  nbytes=*nx*4;
  reclen=nbytes+8;

  if(!(buf = (int *)malloc(reclen))){
    printf("WRTRST: Could not allocate buffer while writing file %s\n",filename);
    exit(1);
  }
 
  if(*zip) {
    (void)ZIPBUFF(filename,rst,    &Open  ,&Form,&Read);
  }
  else {
    if((fd=open(filename,O_RDONLY))<0){
      printf("READRST: Could not open %s error code = %d\n",filename,errno);
      exit(1);
    }
  }


  for (j=0;j<*ny;j++) {

    if(*zip) {
      (void)ZIPBUFF(filename,buf,&reclen,&Form,&Read);
    }
    else {
      if(read(fd,(char *)buf,reclen)<0) {
	printf("READRST: Read error %d on file: %s\n",errno,filename);
	exit(1);
      }
    }

    if(buf[0]!=nbytes) {
      printf("READRST: Wrong record size on file %s: File: %d; Array: %d\n",
	     filename,reclen,*nx);
      exit(1);
    }

    memcpy(rst,(char *)&buf[1],nbytes);
    rst+=*nx;

  }  // Loop over J

  if(*zip) {
    (void)ZIPBUFF(filename,rst,&Close,&Form,&Write);
  }
  else {
    close(fd);
  }

  free(buf);
}

void writerst(int *rst, int *nx, int *ny, char *filename, int *zip)
{
    (void) WRITERST(rst,nx,ny,filename,zip);
}
void writerst_(int *rst, int *nx, int *ny, char *filename, int *zip)
{
    (void) WRITERST(rst,nx,ny,filename,zip);
}
void WRITERST_(int *rst, int *nx, int *ny, char *filename, int *zip)
{
    (void) WRITERST(rst,nx,ny,filename,zip);
}



void readrst(int *rst, int *nx, int *ny, char *filename, int *zip)
{
    (void) READRST(rst,nx,ny,filename,zip);
}
void readrst_(int *rst, int *nx, int *ny, char *filename, int *zip)
{
    (void) READRST(rst,nx,ny,filename,zip);
}
void READRST_(int *rst, int *nx, int *ny, char *filename, int *zip)
{
    (void) READRST(rst,nx,ny,filename,zip);
}

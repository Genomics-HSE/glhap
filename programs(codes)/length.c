// THIS PROGRAM SHOWS LENGTH OF FILES IN LIST

//HOW TO USE
//g++ length.c -lhts -o length
//length
// YOU MUST HAVE A "filelist.txt" FILE IN FOLDER WITH EXECUTIVE FILE WHERE LIST OF FASTA FILES MUST BE WRITTEN MUST BE WRITTEN    

#include <htslib/faidx.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(){
   char *fil;
   fil=(char*)malloc(160);
  void res(char* fain);
  FILE *f = fopen("filelistmt.txt","r");
  while(fscanf(f,"%s",fil)!=EOF)
 res(fil);
  return 0;
}

void res(char* fain){
 faidx_t *fai = NULL;

  fai  = fai_load(fain);
 if(fai==NULL){printf("END");exit(0);}

  assert(faidx_nseq(fai)==1);//file should only contain one entry

  const char *name = faidx_iseq(fai,0);
  int name_len = faidx_seq_len(fai,name);
  printf("%d\n",name_len);
  fai_destroy(fai);
   return;
    }

//TESTREAD
//print letter letter bases of input fastas, then fastas from filelist

//HOW TO USE
//g++ testread.c -lhts -o testread 
//testread f.fa

//ACTUALLY USELESS PROGRAM

#include <htslib/faidx.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(){
   char *fil;
   fil=(char*)malloc(160);
  void res(char* fain);
  FILE *f = fopen("filelist.txt","r");
  while(fscanf(f,"%s",fil)!=EOF)
 res(fil);
  return 0;
}

void res(char* fain){
 faidx_t *fai = NULL;
  printf("input fasta is: %s\n",fain);

  fai  = fai_load(fain);
 if(fai==NULL){printf("END");exit(0);}

  assert(faidx_nseq(fai)==1);//file should only contain one entry

  const char *name = faidx_iseq(fai,0);
  int name_len = faidx_seq_len(fai,name);
  printf("name: %s length: %d\n",name,name_len);

  char *data = fai_fetch(fai,name,&name_len);
  for(int i=0;i<name_len;i++)
    printf("%d\t%c\n",i,data[i]);
    fai_destroy(fai);
   return;
    }

// THIS PROGRAM IS USED FOR CALCULATING NUMBER OF DIFFERENT POSITION BETWEEN  - CONSENSUS FASTA FILE - AND LIST OF CONSENSUS FASTA FILES
//HAMMING LIST COMPARE

// HOW TO USE

//g++ hlistcompare.c -lhts -o hlistcompare
// ./count f1.fa

// YOU MUST HAVE A "filelist.txt" FILE IN FOLDER WITH EXECUTIVE FILE WHERE LIST OF FASTA FILES MUST BE WRITTEN MUST BE WRITTEN 


#include <htslib/faidx.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

int main(int argc,char **argv){
  if(argc==1){
    printf("Supply input fasta\n");
    return 0;
  }
   char *fil;
   fil=(char*)malloc(160);
  void res(char* reference,char* fain);
  FILE *f = fopen("filelist.txt","r");
  char* reference=argv[1];
  while(fscanf(f,"%s",fil)!=EOF)
 res(reference,fil);
  return 0;
}

void res(char* reference,char* fain){
 int j=0;
 faidx_t *fai = NULL;
 faidx_t *ref = NULL;
  fai  = fai_load(fain);
  ref  = fai_load(reference);
  assert(fai!=NULL);
  assert(ref!=NULL);

  assert(faidx_nseq(fai)==1);//file should only contain one entry
  assert(faidx_nseq(ref)==1);//file should only contain one entry

  const char *name = faidx_iseq(fai,0);
  const char *nameref = faidx_iseq(ref,0);
  int name_len = faidx_seq_len(fai,name);
  int nameref_len = faidx_seq_len(ref,nameref);

  char *data = fai_fetch(fai,name,&name_len);
  char *dataref = fai_fetch(ref,nameref,&nameref_len);
  for(int i=0;i<name_len;i++)
    if(toupper(data[i])!=toupper(dataref[i]))j++;
    printf("difference between %s and %s is %d\n",nameref,name,j);
    fai_destroy(fai);
    fai_destroy(ref);
   return;
    }

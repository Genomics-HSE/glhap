//THIS PROGRAM CALCULATES 10 BEST HAPLOGROUP ACCORDING HAMMING DISTANCE
// 
//HOW TO USE 
// 
// g++ hamgrep.c -lhts -o hamgrep
// ./hamgrep f1.fa 
// YOU MUST HAVE A "filelist.txt" FILE IN FOLDER WITH EXECUTIVE FILE WHERE LIST OF FASTA FILES MUST BE WRITTEN MUST BE WRITTEN 
// 
// 

#include <htslib/faidx.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
struct hp
{
  int difference;
  char name[50];
};

int compare(const void * x1,const void * x2)   // функция сравнения элементов массива
{
  return ( ((struct hp*)x1)->difference - ((struct hp*)x2)->difference );              
}

int main(int argc,char **argv){
  if(argc==1){
    printf("Supply input fasta\n");
    return 0;
  }

  printf("haplogroup,distance\n");
  struct hp mt[5438];
   char *fil;
   int i=0;
   fil=(char*)malloc(160);
  void res(char* reference,char* fain, struct hp *mt);
  FILE *f = fopen("filelist.txt","r");
  char* reference=argv[1];
  while(fscanf(f,"%s",fil)!=EOF){
 res(reference,fil,mt+i);i++;}
 qsort(mt,5438,sizeof(mt[0]),compare);
 for(i=0;i<10;i++)
 	printf("%s,%d\n",mt[i].name,mt[i].difference);
 	//printf("%d",compare(mt,mt+2));
  return 0;
}

void res(char* reference,char* fain,struct hp *mt){
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
  //  printf("%d is difference between %s and %s \n",j,nameref,name);
    strcpy(mt->name,name);
    mt->difference=j;
    fai_destroy(fai);
    fai_destroy(ref);
   return;
    }

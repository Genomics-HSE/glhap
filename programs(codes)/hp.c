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
   int tmp,min=99999;
   fil=(char*)malloc(160);
  int dif(char* reference,char* fain);
  FILE *f = fopen("filelist.txt","r");
  char* reference=argv[1];
  while(fscanf(f,"%s",fil)!=EOF){//search minimum
  tmp=dif(reference,fil);
  if(min>tmp)min=tmp;
  }
  rewind(f);//go to the beggining of file
  while(fscanf(f,"%s",fil)!=EOF){
  tmp=dif(reference,fil);
  if(tmp==min)printf("difference between %s and %s is %d\n",reference,fil,tmp);
  }
  return 0;
}

int dif(char* reference,char* fain){
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
    fai_destroy(fai);
    fai_destroy(ref);
   return j;
    }


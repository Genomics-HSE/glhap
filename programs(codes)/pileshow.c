#include <htslib/faidx.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

char * getref(char* str){
    strtok(str,"\t");strtok(NULL,"\t ");return  strtok(NULL,"\t");
}

char * getbase(char* str){
    strtok(NULL,"\t");return  strtok(NULL,"\t");
}

char * getbaq(char* str){
    return  strtok(NULL,"\t");
}

char * getmq(char* str){
    return  strtok(NULL,"\t");
}

int main(int argc,char **argv){
  if(argc==1){
    printf("Supply input mpileup file\n");
    return 0;
  }
   char *fil;
   int num=1,tmp,min=99999;
   fil=(char*)malloc(100000);
   char *ref,*base,*baq,*mq;
  void res(char*,char*,char*,char*);
  FILE *f = fopen(argv[1],"r");
   while(fgets(fil,100000,f)!=NULL){
     printf("chr %d\n",num);
  ref=getref(fil);
  base=getbase(fil);
  baq=getbaq(fil);
  mq=getmq(fil);
  //printf("a is %c",baq[0]);
  res(ref,base,baq,mq);num++;}
  //while(fscanf(f,"%s",fil)!=0)res(3rdword(fil));

puts(baq);
//printf("%lu\n",strlen(base));

  return 0;
}

void res(char* ref,char* base,char* baq, char* mq){
 int j=0,i=0;
for(i;i<strlen(base);i++){
  if(base[i]=='.'||base[i]==',')base[i]=ref[0];
  if(base[i]=='$')continue;
  if(base[i]=='-'){i+=1+base[i+1]-'0';continue;}
if(base[i]=='+'){i+=1+base[i+1]-'0';continue;}
    if(base[i]=='^'){i++;continue;}
      printf("read %10d  base %c \t baq %d \t mq %d\n",j,base[i],baq[j]-33,mq[j]-33);j++;
        }
}

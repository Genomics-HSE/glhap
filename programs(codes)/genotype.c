//THIS PROGRAM CALCULATES GENOTYPE LIKELYHOOD FOR EACH POSITION OF EACH BASE

//HOW TO USE

// g++ likelihood.c -lhts -lz -o likelihood 
// likelihood ./likelihood f.pileup

//THIS PROGRAM HAS OVERFLOW PROBLEM!!!!!!!! 


#include <htslib/faidx.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

long double max(long double x, long double y){
	if(x>y)return x;
	return y;
}

long double pr(int x)
{
	return pow(10,-x/10);
}

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
  puts("chr\ta\tt\tg\tc");
   char *fil;
   int num=1,tmp,min=99999;
   fil=(char*)malloc(100000);
   char *ref,*base,*baq,*mq;
  void res(char*,char*,char*,char*);
  FILE *f = fopen(argv[1],"r");
   while(fgets(fil,100000,f)!=NULL){
     printf("%d\t",num);
  ref=getref(fil);
  base=getbase(fil);
  baq=getbaq(fil);
  mq=getmq(fil);
  //printf("a is %c",baq[0]);
  res(ref,base,baq,mq);num++;}
  //while(fscanf(f,"%s",fil)!=0)res(3rdword(fil));

//puts(baq);
//printf("%lu\n",strlen(base));

  return 0;
}

void res(char* ref,char* base,char* baq, char* mq){
 int j=0,i=0;
 long double a=1,t=1,g=1,c=1;
for(i=0;i<strlen(base);i++){
  if(base[i]=='.'||base[i]==',')base[i]=ref[0];
  if(base[i]=='$')continue;
  if(base[i]=='-'){i+=1+base[i+1]-'0';continue;}
if(base[i]=='+'){i+=1+base[i+1]-'0';continue;}
    if(base[i]=='^'){i++;continue;}
    if(base[i]=='A'||base[i]=='a')
		{
			a*=1-pr(baq[j]-33);
			t*=pr(baq[j]-33)/3;
			g*=pr(baq[j]-33)/3;
			c*=pr(baq[j]-33)/3;
		}
    if(base[i]=='T'||base[i]=='t')
		{
			a*=pr(baq[j]-33)/3;
			t*=1-pr(baq[j]-33);
			g*=pr(baq[j]-33)/3;
			c*=pr(baq[j]-33)/3;
		}
    if(base[i]=='G'||base[i]=='g')
		{
			a*=pr(baq[j]-33)/3;
			t*=pr(baq[j]-33)/3;
			g*=1-pr(baq[j]-33);
			c*=pr(baq[j]-33)/3;
		}
    if(base[i]=='C'||base[i]=='c')
		{
			a*=pr(baq[j]-33)/3;
			t*=pr(baq[j]-33)/3;
			g*=pr(baq[j]-33)/3;
			c*=1-pr(baq[j]-33);
		}

		long double maxi = max(max(a,t),max(g,c));

    j++;

        }

		long double maxi = max(max(a,t),max(g,c));
         printf("%.160Lg\t%.160Lg\t%.160Lg\t%.160Lg\n",logl(a/maxi),logl(t/maxi),logl(g/maxi),logl(c/maxi));
}

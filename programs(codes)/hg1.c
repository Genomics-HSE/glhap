#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <iostream>



void res(char *, char *, char *, char *);

//returns likelihood with of certain base
long double hp(char *, int);

void get_likelihoods(float**,char*);

////////////////////////////////////////////////////////////////////////

long double max(long double x, long double y){
    if(x>y)return x;
    return y;
}

long double atgc[4];//a-0, t-1, g-2, c-3
float hap[16569][4];// all genotype likelihoods for each position and each base
struct haplogroup {
    char name[100];
    long double lh;//sum f likelihoods
};

struct haplogroup mt[5438];

int compare(const void * x1,const void * x2)   // функция сравнения элементов массива
{
    if((((struct haplogroup*)x1)->lh) - (((struct haplogroup*)x2)->lh)>0)return -1;
    if((((struct haplogroup*)x1)->lh) - (((struct haplogroup*)x2)->lh)<0)return 1;
    return 0;
}

float pltoln(int32_t val) {
    return (float)(-val/10);
}

long double pr(int x) {
    return pow(10, -x / 10);
}

char *getref(char *str) {
    strtok(str, "\t");
    strtok(NULL, "\t ");
    return strtok(NULL, "\t");
}

char *getbase(char *str) {
    strtok(NULL, "\t");
    return strtok(NULL, "\t");
}

char *getbaq(char *str) {
    return strtok(NULL, "\t");
}

char *getmq(char *str) {
    return strtok(NULL, "\t");
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
    char fil[100000];
    char *ref, *base, *baq, *mq;
    int i = 0,hnum=0;
    long double sum = 0;
    int num = 0, tmp;


    if (argc == 1) {
        printf("Supply input mpileup file\n");
        return 0;
    }

/////////////////////////////





int *pl = NULL;
    int npl_arr = 0;
    int npl = 0;

    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    int nad_arr = 0;
    int nad     = 0;
    int *ad     = NULL;



    for(int i =0;i<16569;i++)
        for(int j = 0;j<4;j++)
            hap[i][j]=-1;
    if(argc == 2) {
        htsFile *inf = NULL;
        bcf_hdr_t *hdr = NULL;
        bcf1_t *rec = bcf_init();
        inf = bcf_open(argv[1], "r");
        if(inf == NULL) {
            throw std::runtime_error("Unable to open file.");
        }
        hdr = bcf_hdr_read(inf);
        if(hdr == NULL) {
            throw std::runtime_error("Unable to read header.");
        }
        
        std::cout << "haplogroup,log-likelihood" << std::endl;
        while(bcf_read(inf, hdr, rec) == 0) {

            npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);

         //   std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" <<
         //                 rec->pos +1 << "\t";

            bcf_unpack(rec, BCF_UN_ALL);

        
 if (rec->n_allele > 1) {
        for (int i = 0; i < rec->n_allele; ++i)
            for (int j = i; j < rec->n_allele; ++j){
            if (rec->d.allele[i]==rec->d.allele[j]){
            // printf("%s%s %d ",rec->d.allele[i],rec->d.allele[j],pl[i+j]);


            if(strcmp(rec->d.allele[i],"A")==0)hap[rec->pos][0]=-pl[i+j]/10;
            if(strcmp(rec->d.allele[i],"T")==0)hap[rec->pos][1]=-pl[i+j]/10;
            if(strcmp(rec->d.allele[i],"G")==0)hap[rec->pos][2]=-pl[i+j]/10;
            if(strcmp(rec->d.allele[i],"C")==0)hap[rec->pos][3]=-pl[i+j]/10;
            if(strcmp(rec->d.allele[i],"<*>")==0)
                for(int k = 0;k<4;k++)
                    if (hap[rec->pos][k]==-1)
                        hap[rec->pos][k]=-pl[i+j]/10;
                            }
        }}
        }
        bcf_hdr_destroy(hdr);
        bcf_destroy(rec); 
        bcf_close(inf);
    }




///////////////////////////////
   

    FILE *ls = fopen("filelist.txt", "r");
    char haplo[160];
        while (fscanf(ls, "%s", haplo) != EOF) {

        sum = 0;
        faidx_t *fai = NULL;
        fai = fai_load(haplo);
        assert(fai != NULL);
        assert(faidx_nseq(fai) == 1);//file should only contain one entry
        const char *name = faidx_iseq(fai, 0);
        int name_len = faidx_seq_len(fai, name);
        char *data = fai_fetch(fai, name, &name_len);

        for(int j=0;j<16569;j++){
        //printf("data %d %Lg\n", j, hp(data, j));
        sum += hp(data, j);}
       // sum = exp(sum);
      //  printf("log sum of glh for %s is %.30Lf\n", haplo,sum);

        mt[hnum].lh = sum;
        strcpy((mt + hnum) -> name , haplo);
        fai_destroy(fai);
        hnum++;
        
}

//printf("%d\n",hnum);

   qsort(mt,5438 ,sizeof(mt[0]),compare);
    //while(fscanf(f,"%s",fil)!=0)res(3rdword(fil));

    //puts(baq);
    //printf("%lu\n",strlen(base));

//printf("%d",hnum);
    for(i=0;i<10;i++)
        printf("%s,%.30Lg\n",mt[i].name,mt[i].lh);
   /* long double max = 0;
    char * nname;
    for(i=0;i<5438;i++)
        if(max<mt[i].lh){max = mt[i].lh;puts(mt[i].name);}
*/

    return 0;
}



// void res(char *ref, char *base, char *baq, char *mq) {
//     int j = 0, i = 0;
//     long double a = 0, t = 0, g = 0, c = 0;
//     for (i = 0; i < strlen(base); i++) {
//         if (base[i] == '.' || base[i] == ',')base[i] = ref[0];
//         if (base[i] == '$')continue;
//         if (base[i] == '-') {
//             i += 1 + base[i + 1] - '0';
//             continue;
//         }
//         if (base[i] == '+') {
//             i += 1 + base[i + 1] - '0';
//             continue;
//         }
//         if (base[i] == '^') {
//             i++;
//             continue;
//         }
//         if (base[i] == 'A' || base[i] == 'a') {
//             a += logl(1 - pr(baq[j] - 33));
//             t += logl(pr(baq[j] - 33) / 3);
//             g += logl(pr(baq[j] - 33) / 3);
//             c += logl(pr(baq[j] - 33) / 3);
//         }
//         if (base[i] == 'T' || base[i] == 't') {
//             a += logl(pr(baq[j] - 33) / 3);
//             t += logl(1 - pr(baq[j] - 33));
//             g += logl(pr(baq[j] - 33) / 3);
//             c += logl(pr(baq[j] - 33) / 3);
//         }
//         if (base[i] == 'G' || base[i] == 'g') {
//             a += logl(pr(baq[j] - 33) / 3);
//             t += logl(pr(baq[j] - 33) / 3);
//             g += logl(1 - pr(baq[j] - 33));
//             c += logl(pr(baq[j] - 33) / 3);
//         }
//         if (base[i] == 'C' || base[i] == 'c') {
//             a += logl(pr(baq[j] - 33) / 3);
//             t += logl(pr(baq[j] - 33) / 3);
//             g += logl(pr(baq[j] - 33) / 3);
//             c += logl(1 - pr(baq[j] - 33));
//         }

//         j++;

//     }
//     long double maxi = max(max(a,t),max(g,c));

//      atgc[0] = a - maxi;
//     atgc[1] = t - maxi;
//     atgc[2] = g - maxi;
//     atgc[3] = c - maxi;
// }

long double hp(char *hp, int n) {

    switch (hp[n]) {
        case 'a':
            return hap[n][0];
            break;
        case 't':
            return hap[n][1];
            break;
        case 'g':
            return hap[n][2];
            break;
        case 'c':
            return hap[n][3];
            break;
        case 'A':
            return hap[n][0];
            break;
        case 'T':
            return hap[n][1];
            break;
        case 'G':
            return hap[n][2];
            break;
        case 'C':
            return hap[n][3];
            break;
    }
    return 0;
}


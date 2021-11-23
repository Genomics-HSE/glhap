#include <iostream>
#include <stdexcept>
#include <htslib/vcf.h>
#include <string>
#include <cmath>
#include <cstring>


float pltoln(int32_t val) {
    return pow(10.0, -0.1 * val);
}
int usage() {
    std::cerr << "Usage: ./test_htslib_vcfparser example.vcf" << std::endl;
    return 1;
}

int main(int argc, char* argv[]) {

	int *pl = NULL;
 	int npl_arr = 0;
 	int npl = 0;

 	int ngt_arr = 0;
  	int ngt     = 0;
  	int *gt     = NULL;
  	int nad_arr = 0;
  	int nad     = 0;
  	int *ad     = NULL;


  	int pls[16569][4];   //AA TT GG CC;
  	for(int i =0;i<16569;i++)
  		for(int j = 0;j<4;j++)
  			pls[i][j]=-1;
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
        
        std::cout << "chromosome\tposition\tPL" << std::endl;
        while(bcf_read(inf, hdr, rec) == 0) {
//        	 ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);

        	 // std::string c_gt="0/0";
         //  if(gt[0]==2 && gt[1]==4) c_gt="0/1";
         //  if(gt[0]==4 && gt[1]==4) c_gt="1/1";
         //  if(gt[0]==2 && gt[1]==5) c_gt="0|1";
         //  if(gt[0]==4 && gt[1]==5) c_gt="1|1";
         //  if(gt[0]==2 && gt[1]==3) c_gt="0|0";
         //  if(gt[0]==4 && gt[1]==3) c_gt="1|0";

        	npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);

            std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" <<
                          rec->pos +1 << "\t";
                          //std::couut << npl_arr
            // for(int i =0;i<npl;i++)
            	// std::cout << pl[i]<< " \t ";



            bcf_unpack(rec, BCF_UN_ALL);

		
 if (rec->n_allele > 1) {
        for (int i = 0; i < rec->n_allele; ++i)
        	for (int j = i; j < rec->n_allele; ++j){
          	if (rec->d.allele[i]==rec->d.allele[j]){
             printf("%s%s %d ",rec->d.allele[i],rec->d.allele[j],pl[i+j]);


            if(strcmp(rec->d.allele[i],"A")==0)pls[rec->pos][0]=pl[i+j];
            if(strcmp(rec->d.allele[i],"T")==0)pls[rec->pos][1]=pl[i+j];
            if(strcmp(rec->d.allele[i],"G")==0)pls[rec->pos][2]=pl[i+j];
            if(strcmp(rec->d.allele[i],"C")==0)pls[rec->pos][3]=pl[i+j];
            if(strcmp(rec->d.allele[i],"<*>")==0)
            	for(int k = 0;k<4;k++)
  					if (pls[rec->pos][k]==-1)
  						pls[rec->pos][k]=pl[i+j];
  							}
        }}
        for(int k = 0;k<4;k++)
        	printf("%d  ",pls[rec->pos][k]);
        printf("\n");	


            	std::cout<< "\n";
           // std::cout << "\t\tGL\t";
            //for(int i =0;i<npl_arr;i++)
            //	std::cout << pltoln(pl[i])<< "\t";
            std::cout << "\n";
        //                  rec->qual << "\t" <<
        //                  std::endl;
        // }

        }
        bcf_hdr_destroy(hdr);
        bcf_destroy(rec); 
        bcf_close(inf);
    } else {
        return usage();
    }



    // for(int i =0;i<16569;i++){
  		// for(int j = 0;j<4;j++)
  		// 	printf("%d ",pls[i][j]);
  		// printf("\n");
  	// }
    return 0;
}

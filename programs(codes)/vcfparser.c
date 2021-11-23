#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

void main(int argc,char **argv)
{
    char* fname = argv[1];

    //open vcf file
    htsFile *fp    = hts_open(fname,"rb");
    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);

    bcf1_t *rec    = bcf_init();

    //save for each vcf record
    while ( bcf_read(fp, hdr, rec)>=0 )
    {
        //unpack
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);


        //read CHROM
        printf("%s\n", bcf_hdr_id2name(hdr, rec->rid));
        //readqual
        printf("%s\n", bcf_hdr_id2int(hdr, BCF_DT_ID, "PL"));
    }



    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
}
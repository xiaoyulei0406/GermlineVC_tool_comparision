import "variant_calling.wdl" as vartool
import "germline_filter.wdl" as filtertool
workflow hardfilter {
	
	File input_bam
	String output_prefix
	File knowsite1_bqsr_golden_indel
	File knowsite2_bqsr_phase1_SNPs
	File knowsite3_bqsr_dbsnp
	File dbsnp
	File Ref
	File bedFile


	call vartool.bqsr as bqsr {
		input:
			bam = input_bam,
			Ref_file = Ref,
			knowsite1 = knowsite1_bqsr_golden_indel,
			knowsite2 = knowsite2_bqsr_phase1_SNPs,
			knowsite3 = knowsite3_bqsr_dbsnp,
			sample_name_flag = output_prefix
	}
	call vartool.haplotypecaller {
		input:
			bam = bqsr.recalibrated_bam,
			sample_name_flag = output_prefix,
			Ref_file = Ref, dbSNP = dbsnp,
			recal_table = bqsr.recal_table,
			intervals = bedFile
	}
	call filtertool.vf_snp as vf_snp {
		input: Ref_file = Ref, vcf = haplotypecaller.vcf
	}
	call filtertool.vf_indel as vf_indel {
		input: Ref_file = Ref, vcf = haplotypecaller.vcf
	}
	call filtertool.merge_snp_indel {
		input: Ref_file = Ref, indels = vf_indel.indels, snps = vf_snp.snps
	}
    output {
        File snp = vf_snp.snps
        File indel = vf_indel.indels
        File combine_vcf = merge_snp_indel.filter
    }
}

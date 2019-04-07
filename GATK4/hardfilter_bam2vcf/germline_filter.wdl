task vf_snp {
	File Ref_file
	File vcf
#	String outdir

	command {
		java -jar /gatk/gatk.jar SelectVariants -R ${Ref_file} -V ${vcf} \
		--select-type-to-include SNP -O raw_snp.vcf

		java -jar /gatk/gatk.jar VariantFiltration -R ${Ref_file} -V raw_snp.vcf \
		--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
		--filter-name "HQ_SNP_filter" -O germline_snp.vcf

		rm raw_snp.vcf
		rm raw_snp.vcf.idx
	}
	output { 
		File snps = "germline_snp.vcf"
	}

	runtime {
		docker:"gatk4/wd:latest"
		cpu:"4"
		memory:"30G"
	}
}

task vf_indel {
	File Ref_file
	File vcf
#	String outdir

	command {
		java -jar /gatk/gatk.jar SelectVariants -R ${Ref_file} -V ${vcf} \
		--select-type-to-exclude SNP -O raw_indel.vcf

		java -jar /gatk/gatk.jar VariantFiltration -R ${Ref_file} -V raw_indel.vcf \
		--filter-expression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0" \
		--filter-name "HQ_Indel_filter" -O germline_indel.vcf

		rm raw_indel.vcf
		rm raw_indel.vcf.idx
	}
	output { 
		File indels = "germline_indel.vcf"
	}

	runtime {
		docker:"gatk4/wd:latest"
		cpu:"4"
		memory:"30G"
	}
}


task merge_snp_indel {
	File Ref_file
	File indels
	File snps
#	String outdir
	command {
		vcfcombine ${snps} ${indels} > combine.vcf
	}
	output { 
		File filter = "combine.vcf"
	}
	runtime {
		docker:"erictdawson/svdocker:latest"
		cpu:"4"
		memory:"30G"
	}

	meta {
		author: "Chunlei Yu"
	}
}

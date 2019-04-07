task bqsr {
	
	File input_bam
	File reference_fasta
	File knowsite1
	File knowsite2
	File knowsite3
	String sample_name_flag

	command <<<
		gatk BaseRecalibrator --input ${input_bam} --output ${sample_name_flag}.bqsrtable --reference ${reference_fasta} --known-sites ${knowsite1} --known-sites ${knowsite2} --known-sites ${knowsite3}

		gatk ApplyBQSR -R ${reference_fasta} -I ${input_bam} -bqsr ${sample_name_flag}.bqsrtable -O ${sample_name_flag}.bqsr.bam
	>>>

	output {
		File recal_table = "${sample_name_flag}.bqsrtable"
		File recalibrated_bam = "${sample_name_flag}.bqsr.bam"
		File outbai = "${sample_name_flag}.bqsr.bam.bai"
	}
	
	runtime{
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"50G"
	}
}

task RunHC4 {
	File input_bam
	File reference_fasta
	String output_prefix
	File interval_list
	String extra_args

	command {
		gatk HaplotypeCaller \
		-R ${reference_fasta} \
		-I ${input_bam} \
		-O ${output_prefix}_hc4.vcf.gz \
		-L ${interval_list} \
		-bamout ${output_prefix}_bamout.bam \
		${extra_args}
	}

	output {
		File bamout = "${output_prefix}_bamout.bam"
		File bamout_index = "${output_prefix}_bamout.bai"
		File raw_vcf = "${output_prefix}_hc4.vcf.gz"
		File raw_vcf_index = "${output_prefix}_hc4.vcf.gz.tbi"
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"10"
		memory:"30G"
	}
}

task GenomicsDBImport {
	Array[File] GVCFs
	String flag
#	String outdir
	File bedfile

	command <<<
		gatk GenomicsDBImport \
		-V ${sep=' -V ' GVCFs} \
		--genomicsdb-workspace-path ${flag} \
    	--intervals ${bedfile}

		tar -cf ${flag}.tar ${flag}
	>>>
	output { 
		File db = "${flag}.tar"
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"30G"
	}
}

task GenotypeGVCFs {
	File Ref_file
#	Array[File] GVCFs
	File my_database
#	String outdir
	File dbSNP
	String sample_name_flag

	command <<<
		tar -xf ${my_database}
		WORKSPACE=$( basename ${my_database} .tar)
		echo $WORKSPACE

		gatk GenotypeGVCFs -R ${Ref_file} \
		-V gendb://$WORKSPACE \
		-G StandardAnnotation \
		-D ${dbSNP} \
		-O ${sample_name_flag}.vcf
	>>>

	output { 
		File vcf = "${sample_name_flag}.vcf"
		File output_vcf_index = "${sample_name_flag}.tbi"
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"30G"
	}
}

task MakeSitesOnlyVcf {
	File input_vcf
	String sample_name_flag

	command {
		java -jar /opt/picard/picard.jar \
		MakeSitesOnlyVcf \
		INPUT=${input_vcf} \
		OUTPUT=${sample_name_flag}.sites_only.vcf
	}
	runtime {
		docker:"aws_tool:2.1.4"
		cpu:"4"
		memory:"30G"
	}
	output { 
		File sites_only_vcf  = "${sample_name_flag}.sites_only.vcf"
		File sites_only_vcf_index = "${sample_name_flag}.sites_only.vcf.tbi"
	}
}

task VariantRecalibratorSNP {
	File input_vcf
	File knowsite2_phase1_SNPs
	File knowsite3_hapmap
	File knowsite4_omni
	File knowsite6_bqsr_dbsnp
	File reference_fasta
	String output_prefix

	command{
		gatk VariantRecalibrator \
		-R ${reference_fasta} \
		-V ${input_vcf} \
		--mode SNP \
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${knowsite3_hapmap} \
		-resource:omni,known=false,training=true,truth=true,prior=12.0 ${knowsite4_omni} \
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${knowsite2_phase1_SNPs} \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${knowsite6_bqsr_dbsnp} \
		-an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
		-tranche 100.0 -tranche 99.9 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		--tranches-file ${output_prefix}.tranches \
		-O ${output_prefix}.recal
	}
	output {
		File SNP_recal = "${output_prefix}.recal"
		File SNP_recalIndex = "${output_prefix}.recal.idx"
		File SNP_tranches = "${output_prefix}.tranches"
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"3"
		memory:"30G"
	}
}
task ApplyRecalibrationSNP {
	File reference_fasta
	File input_vcf
	File SNP_tranches
	File SNP_recal
	String output_prefix

	command {
		gatk ApplyVQSR \
		-V ${input_vcf} \
		-R ${reference_fasta} \
		--mode SNP \
		-ts-filter-level 99.0 \
		-recal-file ${SNP_recal} \
		-tranches-file ${SNP_tranches} \
		-O ${output_prefix}.SNP.g.vcf

	}
	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"3"
		memory:"30G"
	}
	output {
		File snp_vcf = "${output_prefix}.SNP.g.vcf"
		File snp_vcf_index = "${output_prefix}.SNP.g.vcf.idx"
	}
}
task VariantRecalibratorINDEL {
	File input_vcf
	File knowsite1_golden_indel
	File knowsite6_bqsr_dbsnp
	File reference_fasta
	String output_prefix

	command {
		gatk VariantRecalibrator \
		-R ${reference_fasta} \
		-V ${input_vcf} \
		--mode INDEL \
		-resource:mills,known=false,training=true,truth=true,prior=12.0 ${knowsite1_golden_indel}  \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${knowsite6_bqsr_dbsnp} \
		-an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		--max-gaussians 4 \
		--tranches-file ${output_prefix}.tranches \
		-O ${output_prefix}.recal
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"3"
		memory:"30G"
	}
	output {
		File indel_recal = "${output_prefix}.recal"
		File indel_recalIndex = "${output_prefix}.recal.idx"
		File indel_tranches = "${output_prefix}.tranches"
	}
}

task ApplyRecalibrationINDEL {
	File reference_fasta
	File input_vcf
	File indel_tranches
	File indel_recal
	String output_prefix

	command {
		gatk ApplyVQSR \
		-V ${input_vcf} \
		-R ${reference_fasta} \
		--mode INDEL \
		-ts-filter-level 95.0 \
		--tranches-file ${indel_tranches} \
		--recal-file ${indel_recal} \
		-O ${output_prefix}.indel.g.vcf
	}
	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"30G"
	}
	output {
		File indel_vcf = "${output_prefix}.indel.g.vcf"
		File indel_vcf_index = "${output_prefix}.indel.g.vcf.idx"
    }
}


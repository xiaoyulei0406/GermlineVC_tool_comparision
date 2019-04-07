task bqsr {
	
	File input_bam
	File reference_fasta
	File knowsite1
	File knowsite2
	File knowsite3
	String sample_name_flag

	command <<<
		java -jar /gatk/gatk.jar BaseRecalibrator --input ${input_bam} --output ${sample_name_flag}.bqsrtable --reference ${reference_fasta} --known-sites ${knowsite1} --known-sites ${knowsite2} --known-sites ${knowsite3}

		java -jar /gatk/gatk.jar ApplyBQSR -R ${reference_fasta} -I ${input_bam} -bqsr ${sample_name_flag}.bqsrtable -O ${sample_name_flag}.bqsr.bam
	>>>

	output {
		File recal_table = "${sample_name_flag}.bqsrtable"
		File recalibrated_bam = "${sample_name_flag}.bqsr.bam"
		File outbai = "${sample_name_flag}.bqsr.bam.bai"
	}
	
	runtime{
		docker:"gatk4/wd:latest"
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

task CNNScoreVariants {
	File input_vcf
	File input_bam
	File reference_fasta
	String output_prefix
	File interval_list
    String tensor_type

	command{
		gatk CNNScoreVariants \
		-I ${input_bam} \
		-R ${reference_fasta} \
		-V ${input_vcf} \
		-O ${output_prefix}_cnn_annotated.vcf.gz \
		-L ${interval_list} \
		${"--tensor-type " + tensor_type}
	}
	output {
#		Array[File] log = glob("gatkStreamingProcessJournal*")
		File cnn_annotated_vcf = "${output_prefix}_cnn_annotated.vcf.gz"
		File cnn_annotated_vcf_index = "${output_prefix}_cnn_annotated.vcf.gz.tbi"
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"10"
		memory:"30G"
	}
}

task MergeVCFs {
	Array[File] input_vcfs
	String output_prefix

	command {
		gatk MergeVcfs \
		-I ${sep=' -I ' input_vcfs} \
		-O ${output_prefix}._cnn_scored.vcf.gz
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"30G"
	}
	output {
        File merged_vcf = "${output_prefix}._cnn_scored.vcf.gz"
        File merged_vcf_index = "${output_prefix}._cnn_scored.vcf.gz.tbi"
    }
}

task FilterVariantTranches {
	File input_vcf
	File knowsite1_golden_indel
    File knowsite2_phase1_SNPs
    File knowsite3_hapmap
    File knowsite4_omni
    File knowsite5_phase1_indels 
    File knowsite6_bqsr_dbsnp
	String output_prefix
	String snp_tranches
	String indel_tranches
	String info_key

	command {
		gatk FilterVariantTranches \
		-V ${input_vcf} \
		--output ${output_prefix}_cnn_filtered.vcf.gz \
		-resource ${knowsite1_golden_indel} \
		-resource ${knowsite2_phase1_SNPs} \
		-resource ${knowsite3_hapmap} \
		-resource ${knowsite4_omni} \
		-resource ${knowsite5_phase1_indels} \
		-resource ${knowsite6_bqsr_dbsnp} \
		--info-key ${info_key} \
		${snp_tranches} ${indel_tranches}
	}
	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"30G"
	}
	output {
		File cnn_filtered_vcf = "${output_prefix}_cnn_filtered.vcf.gz"
		File cnn_filtered_vcf_index ="${output_prefix}_cnn_filtered.vcf.gz.tbi"
	}
}

task SamtoolsMergeBAMs {
	Array[File] input_bams
	String output_prefix

	command {
		samtools merge ${output_prefix}_bamout.bam ${sep=' ' input_bams}
		samtools index ${output_prefix}_bamout.bam ${output_prefix}_bamout.bai
	}

	runtime {
		docker:"broadinstitute/gatk:latest"
		cpu:"4"
		memory:"30G"
	}

	output {
		File bamout = "${output_prefix}_bamout.bam"
		File bamout_index = "${output_prefix}_bamout.bai"
	}
}

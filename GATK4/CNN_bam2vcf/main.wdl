# This workflow takes an input BAM to call variants with HaplotypeCaller
# Then filters the calls with the CNNVariant neural net tool
# The site-level scores are added to the INFO field of the VCF.
# The architecture arguments, info_key and tensor type arguments MUST be in agreement
# (e.g. 2D models must have tensor_type of read_tensor and info_key CNN_2D, 1D models have tensor_type reference and info_key CNN_1D)
# The INFO field key will be "1D_CNN" or "2D_CNN" depending on the neural net architecture used for inference.
# The architecture arguments specify pre-trained networks.
# New networks can be trained by the GATK tools: CNNVariantWriteTensors and CNNVariantTrain
# The bam could be generated by the single-sample pipeline
# (https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl)
# Also accepts a BAM as the input file in which case a BAM index is required as well.

import "cnn_variant_common_tasks.wdl" as CNNTasks

workflow CNN {
    File input_bam                  # Aligned BAM files
    File? input_bam_index           # Index for an aligned BAM file
    File reference_fasta
    String output_prefix             # Identifying string for this run will be used to name all output files
    String? tensor_type              # What kind of tensors the Neural Net expects (e.g. reference, read_tensor)
    String info_key                  # The score key for the info field of the vcf (e.g. CNN_1D, CNN_2D)
    String snp_tranches              # Filtering threshold(s) for SNPs in terms of sensitivity to overlapping known variants in resources
    String indel_tranches            # Filtering threshold(s) for INDELs in terms of sensitivity to overlapping known variants in resources
    File calling_intervals
    Int scatter_count                # Number of shards for parallelization of HaplotypeCaller and CNNScoreVariants
    String extra_args                # Extra arguments for HaplotypeCaller
    String output_prefix

    File knowsite1_golden_indel      # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
    File knowsite2_phase1_SNPs       # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
    File knowsite3_hapmap            # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
    File knowsite4_omni             # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
    File knowsite5_phase1_indels    # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)
    File knowsite6_bqsr_dbsnp       # File of VCF file names of resources of known SNPs and INDELs, (e.g. mills, gnomAD)

    call CNNTasks.bqsr {
        input:
            input_bam = input_bam,
            reference_fasta = reference_fasta,
            knowsite1 = knowsite1_golden_indel,
            knowsite2 = knowsite2_phase1_SNPs,
            knowsite3 = knowsite6_bqsr_dbsnp,
            sample_name_flag = output_prefix,


    }
    call CNNTasks.RunHC4 {
        input:
            input_bam = bqsr.recalibrated_bam,
            reference_fasta = reference_fasta,
            output_prefix = output_prefix,
            interval_list = calling_intervals,
            extra_args = extra_args
    }

     call CNNTasks.CNNScoreVariants {
        input:
            input_vcf = RunHC4.raw_vcf,
	    input_bam = RunHC4.bamout,
            reference_fasta = reference_fasta,
            tensor_type = tensor_type,        
            output_prefix = output_prefix,
            interval_list = calling_intervals
    }


    call CNNTasks.FilterVariantTranches {
        input:
            input_vcf = CNNScoreVariants.cnn_annotated_vcf,
            output_prefix = output_prefix,
            snp_tranches = snp_tranches,
            indel_tranches = indel_tranches,
            info_key = info_key,
            knowsite1_golden_indel = knowsite1_golden_indel,
            knowsite2_phase1_SNPs = knowsite2_phase1_SNPs,
            knowsite3_hapmap = knowsite3_hapmap,
            knowsite4_omni = knowsite4_omni,
            knowsite5_phase1_indels = knowsite5_phase1_indels,
            knowsite6_bqsr_dbsnp = knowsite6_bqsr_dbsnp
    }

    output {
        File cnn_filtered_vcf = FilterVariantTranches.cnn_filtered_vcf
        File cnn_filtered_vcf_index =FilterVariantTranches.cnn_filtered_vcf_index
    }

}

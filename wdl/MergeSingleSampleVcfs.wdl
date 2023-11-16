###############################
#    RAS Modifiers Project    #
###############################

# Copyright (c) 2023-Present Ryan L. Collins and the Van Allen/Gusev/Haigis Laboratories
# Distributed under terms of the GNU GPL v2.0 License (see LICENSE)
# Contact: Ryan L. Collins <Ryan_Collins@dfci.harvard.edu>

# Merge multiple single-sample VCFs into a joint VCF using bcftools


version 1.0


workflow MergeSingleSampleVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String prefix
    String docker = "vanallenlab/g2c_pipeline:latest"
  }

    call MergeVcfs {
      input:
        vcfs = vcfs,
        vcf_idxs = vcf_idxs,
        prefix = prefix,
        docker = docker
  }
  
  output {
    File merged_vcf = MergeVcfs.merged_vcf
    File merged_vcf_idx = MergeVcfs.merged_vcf_idx
  }
}


task MergeVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String prefix

    Int n_cpu = 4
    Float mem_gb_per_cpu = 1.75
    Int? disk_gb

    String docker = "vanallenlab/g2c_pipeline:latest"
  }

  Int default_disk_gb = ceil(3 * size(vcfs, "GB")) + 25
  Int use_disk_gb = select_first([disk_gb, default_disk_gb])
  Float mem_gb = n_cpu * mem_gb_per_cpu

  command <<<
    set -eu -o pipefail

    bcftools merge \
      --missing-to-ref \
      --no-version \
      --threads ~{n_cpu} \
      -Oz -o ~{prefix}.vcf.gz \
      --file-list ~{write_lines(vcfs)}

    tabix -f -p vcf ~{prefix}.vcf.gz
  >>>

  output {
    File merged_vcf = "~{prefix}.vcf.gz"
    File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
  }

  runtime {
    cpu: n_cpu
    memory: mem_gb + " GiB"
    disks: "local-disk " + use_disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}
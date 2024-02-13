//
// Optional sub-workflow for extracting and checking fingerprints
//
include { PICARD_CREATESEQUENCEDICTIONARY           } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_EXTRACTFINGERPRINT                 } from '../../../modules/nf-core/picard/extractfingerprint/main'
include { PICARD_CROSSCHECKFINGERPRINTS             } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'

workflow FINGERPRINTS_PICARD {

    take:
    bam
    fasta
    haplotype_map
    crosscheckfingerprint_sample_map
    fasta_fai

    main:

    ch_versions = Channel.empty()

    // Create dictionary

    PICARD_CREATESEQUENCEDICTIONARY ( fasta )
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions.first())

    // Extract fingerprint from Picard markduplicates output

    PICARD_EXTRACTFINGERPRINT ( bam, haplotype_map, fasta, fasta_fai, PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict )
    ch_versions = ch_versions.mix(PICARD_EXTRACTFINGERPRINT.out.versions.first())

    ch_merged_vcfs = PICARD_EXTRACTFINGERPRINT
    .out
    .vcf
    .map { it ->
        it[1]
    }
    .collect()
    .map { it -> [[id:'fingerprints'], it]
    }

    // Cross-check fingerprints

    ch_input_2 = Channel.fromPath("${projectDir}/assets/dummy_file.txt")

    PICARD_CROSSCHECKFINGERPRINTS ( ch_merged_vcfs, ch_input_2, haplotype_map, crosscheckfingerprint_sample_map )
    ch_versions = ch_versions.mix(PICARD_CROSSCHECKFINGERPRINTS.out.versions.first())

    emit:
    vcf                     = PICARD_EXTRACTFINGERPRINT.out.vcf           // channel: [ val(meta), [ vfc ] ]
    crosscheck_metrics      = PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics          // channel: [ val(meta), [ crosscheck_metrics ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}


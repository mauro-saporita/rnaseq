//
// Optional sub-workflow for extracting and checking fingerprints
//
include { PICARD_CREATESEQUENCEDICTIONARY           } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { PICARD_EXTRACTFINGERPRINT                 } from '../../../modules/local/picard_extractfingerprint/main'
include { PICARD_CROSSCHECKFINGERPRINTS             } from '../../../modules/nf-core/picard/crosscheckfingerprints/main'

workflow FINGERPRINTS_PICARD {

    take:
    fasta
    haplotype_map
    fasta_fai

    main:

    ch_versions = Channel.empty()

    // Create dictionary

    PICARD_CREATESEQUENCEDICTIONARY ( fasta )
    ch_versions = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions.first())

    // Extract fingerprint from Picard markduplicates output

    PICARD_EXTRACTFINGERPRINT ( PICARD_MARKDUPLICATES.out.bam, haplotype_map, fasta, fasta_fai, PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict )
    ch_versions = ch_versions.mix(PICARD_EXTRACTFINGERPRINT.out.versions.first())

    // Cross-check fingerprints

    PICARD_CROSSCHECKFINGERPRINTS ( PICARD_EXTRACTFINGERPRINT.out.vcf, haplotype_map )
    ch_versions = ch_versions.mix(PICARD_CROSSCHECKFINGERPRINTS.out.versions.first())

    emit:
    vcf                     = PICARD_EXTRACTFINGERPRINT.out.vcf           // channel: [ val(meta), [ vfc ] ]
    crosscheck_metrics      = PICARD_CROSSCHECKFINGERPRINTS.out.crosscheck_metrics          // channel: [ val(meta), [ crosscheck_metrics ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}


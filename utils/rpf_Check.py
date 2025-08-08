#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Check.py


from utils.ribo import ArgsParser
from utils.ribo.Quality import *


def main():
    ArgsParser.now_time()
    print('\nCheck the RPFs mapping condition.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_bam_check_parser()
    rpf_quality = Quality(args)

    print('\nStep2: Import the transcripts annotation.', flush=True)
    rpf_quality.read_transcript()

    print('\nStep3: Import the bam file.', flush=True)
    rpf_quality.sort_index_bam()
    rpf_quality.fliter_mrna_reads()

    print('\nStep4: Sort and index the bam file.', flush=True)
    rpf_quality.merge_sort_index_bam()

    print('\nStep5: Detect the type of sequence profile.', flush=True)
    if not rpf_quality.profile:
        rpf_quality.detect_seq_type()
    else:
        print("{bam} is specified as {ribo}-seq.".format(bam=rpf_quality.sample_file, ribo=rpf_quality.profile),
              flush=True)

    print('\nStep6: Summary the length distribution of reads aligned to mRNA.', flush=True)
    rpf_quality.write_length_distr()

    if args.saturation:
        print('\nStep7: Check the RPFs saturation.', flush=True)
        rpf_quality.rpf_saturation()
        rpf_quality.draw_gene_saturation()
        rpf_quality.draw_rpf_saturation()
    else:
        print('\nStep7: Don not check the RPFs saturation.', flush=True)

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

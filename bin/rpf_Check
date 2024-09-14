#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Check.py


from foo import ArgsParser
from foo.Quality import *


def main():
    ArgsParser.now_time()
    print('Check the RPFs mapping condition.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.rpf_bam_check_parser()
    rpf_quality = Quality(args)

    print('Step2: Import the transcripts annotation.\n', flush=True)
    rpf_quality.read_transcript()

    print('Step3: Import the bam file.\n', flush=True)
    rpf_quality.sort_index_bam()
    rpf_quality.fliter_mrna_reads()

    print('Step4: Sort and index the bam file.\n', flush=True)
    rpf_quality.merge_sort_index_bam()

    print('Step5: Detect the type of sequence profile.\n', flush=True)
    if not rpf_quality.profile:
        rpf_quality.detect_seq_type()
    else:
        print("{bam} is specified as {ribo}-seq.\n".format(bam=rpf_quality.sample_file, ribo=rpf_quality.profile),
              flush=True)

    print('Step6: Summary the length distribution of reads aligned to mRNA.\n', flush=True)
    rpf_quality.write_length_distr()

    if args.saturation:
        print('Step7: Check the RPFs saturation.\n', flush=True)
        rpf_quality.rpf_saturation()
        rpf_quality.draw_gene_saturation()
        rpf_quality.draw_rpf_saturation()
    else:
        print('Step7: Don not check the RPFs saturation.\n', flush=True)

    print('All done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

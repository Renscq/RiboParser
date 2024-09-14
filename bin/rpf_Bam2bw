#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboparser
# @Script  : rpf_bam2bw.py


from foo import ArgsParser


def main():
    ArgsParser.now_time()
    print('\nConvert genome bam reads to bedgraph.\n', flush=True)
    print('Step1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_bam2bw_args_parser()

    from foo import Bam2Wig
    bam_attr = Bam2Wig.Bam2Wig(args)

    print('\nStep2: Import the P-site offset.', flush=True)
    bam_attr.read_offset()

    print('\nStep3: Import the alignment bam file.', flush=True)
    bam_attr.set_tag_num()
    bam_attr.import_bam()

    print('\nStep4: convert bam reads to density.', flush=True)
    bam_attr.convert_dict_to_dataframe()

    print('\nStep5: Normalize the reads density.', flush=True)
    bam_attr.norm_rpm()

    print('\nStep6: Output reads density.', flush=True)
    bam_attr.output_bed()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

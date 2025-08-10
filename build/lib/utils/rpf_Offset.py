#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : detect_offset.py


from utils.ribo import ArgsParser
from utils.ribo.Offset import *


def main():
    ArgsParser.now_time()
    print('\nDetect the p-site offset.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.offset_args_parser()
    offset_attr = Offset(args)

    print('\nStep2: Import the transcripts annotation.', flush=True)
    offset_attr.read_transcript()

    print('\nStep3: Import the bam file.', flush=True)
    offset_attr.get_mrna_reads()

    # if args.mode == 'SSCBM':
    print('\nStep4: Detect the SSCBM offset of sequence profile.', flush=True)
    offset_attr.get_tis_offset()
    offset_attr.adjust_tis_offset()
    offset_attr.write_tis_offset()

    offset_attr.draw_tis_heatmap()

    print('\nStep5: Detect the RSBM offset of sequence profile.', flush=True)
    offset_attr.get_frame_offset()
    offset_attr.format_frame_offset()
    offset_attr.adjust_frame_offset()
    offset_attr.write_frame_offset()

    offset_attr.draw_frame_heatmap()

    # elif args.mode == 'RSBM':
    #     print('Step4: Detect the RSBM offset of sequence profile.\n', flush=True)
    #     offset_attr.get_frame_offset()
    #     offset_attr.format_frame_offset()
    #     offset_attr.adjust_frame_offset()
    #     offset_attr.write_frame_offset()

    #     offset_attr.draw_frame_heatmap()

    ArgsParser.now_time()
    print('\nAll done.', flush=True)


if __name__ == '__main__':
    main()

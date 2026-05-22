#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Offset_RSBM.py


from utils.ribo import ArgsParser
from utils.ribo.Offset_RSBM import *


def main():
    ArgsParser.now_time()
    print('\nDetect the p-site offset.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.offset_rsbm_args_parser()
    offset_attr = Offset(args)

    print('\nStep2: Import the transcripts annotation.', flush=True)
    offset_attr.read_transcript()

    print('\nStep3: Import the bam file.', flush=True)
    offset_attr.get_mrna_reads()

    print('\nStep4: Detect the RSBM offset of sequence profile.', flush=True)
    offset_attr.get_frame_offset()
    offset_attr.format_frame_offset()
    offset_attr.adjust_frame_offset()

    print('\nStep5: Output the frame offset.', flush=True)
    offset_attr.write_frame_offset()

    print('\nStep6: Draw the frame offset heatmap.', flush=True)
    offset_attr.draw_frame_heatmap()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

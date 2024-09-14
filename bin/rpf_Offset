#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : detect_offset.py


from foo import ArgsParser
from foo.Offset_mono import *


def main():
    ArgsParser.now_time()
    print('Detect the p-site offset.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.offset_args_parser()
    offset_attr = Offset(args)

    print('Step2: Import the transcripts annotation.\n', flush=True)
    offset_attr.read_transcript()

    print('Step3: Import the bam file.\n', flush=True)
    offset_attr.get_mrna_reads()

    print('Step4: Detect the tis offset of sequence profile.\n', flush=True)

    offset_attr.get_tis_offset()
    offset_attr.adjust_tis_offset()
    offset_attr.write_tis_offset()
    offset_attr.draw_tis_heatmap()

    print('Step5: Detect the frame offset of sequence profile.\n', flush=True)
    offset_attr.get_frame_offset()
    offset_attr.format_frame_offset()
    offset_attr.adjust_frame_offset()

    offset_attr.write_frame_offset()
    offset_attr.draw_frame_heatmap()

    ArgsParser.now_time()
    print('All done.\n', flush=True)


if __name__ == '__main__':
    main()

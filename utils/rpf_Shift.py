#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Shift.py


from utils.ribo import Shift

import argparse
from utils.ribo.ArgsParser import args_print, file_check, now_time


def frame_shift_args_parser():
    parser = argparse.ArgumentParser(description="This script is used to draw the frame shift plot.")

    # arguments for the Required arguments
    input_group = parser.add_argument_group('Required arguments')
    input_group.add_argument('-r', dest="rpf", required=True, type=str,
                             help="the name of input RPFs file in TXT format.")
    input_group.add_argument(
        '-o', dest="output", required=True, type=str,
        help="the prefix of output file."
    )

    # arguments for the ribo-seq parsing
    input_group.add_argument('-t', dest="transcript", required=False, type=str,
                             help="the name of input transcript filein TXT format.")
    parser.add_argument('-p', dest="period", required=False, type=int, default=45,
                        help="the minimum in-frame value for frame shifting screen, range [0 - 1]. (default: %(default)s).")
    parser.add_argument('-m', dest="min", required=False, type=int, default=50,
                        help="retain transcript with more than minimum RPFs. (default: %(default)s).")
    parser.add_argument('--tis', dest="tis", required=False, type=int, default=0,
                        help="the number of codons after TIS will be discarded.. (default: %(default)s AA).")
    parser.add_argument('--tts', dest="tts", required=False, type=int, default=0,
                        help="the number of codons before TTS will be discarded.. (default: %(default)s AA).")
    args = parser.parse_args()
    file_check(args.rpf)
    args_print(args)

    return args

def main():
    now_time()
    print('\nDraw the frame shifting plot.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = frame_shift_args_parser()

    print('\nStep2: Import the RPFs file.', flush=True)
    rpfs = Shift.Shift(args)
    rpfs.import_rpf()

    print('\nStep3: Calculate the 3nt periodicity.', flush=True)
    rpfs.calc_3nt_period()

    print('\nStep4: Ouput the 3nt periodicity.', flush=True)
    rpfs.output_meta()

    print('\nStep5: Filter frame shift.', flush=True)
    rpfs.filter_frame_shift()

    print('\nStep6: Output the frame shift.', flush=True)
    rpfs.output_frame_shift()

    print('\nStep7: Draw the frame shifting plot.', flush=True)
    rpfs.draw_frame_shift_count()

    print('\nAll done.\n', flush=True)
    now_time()


if __name__ == '__main__':
    main()

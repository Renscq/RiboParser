#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_Shift.py


from utils.ribo import ArgsParser
from utils.ribo import Shift

def main():
    ArgsParser.now_time()
    print('\nDraw the frame shifting plot.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.frame_shift_args_parser()

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
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

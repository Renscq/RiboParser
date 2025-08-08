#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_pausing.py


from utils.ribo import ArgsParser
from utils.ribo import Pausing


def main():
    ArgsParser.now_time()
    print('\nCalculate the relative codon pausing score.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.rpf_pausing_args_parser()
    pausing_score = Pausing.Pausing(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    pausing_score.read_rpf()

    print('\nStep3: Calculate the ribosome pausing score.', flush=True)
    pausing_score.get_pausing_score()

    print('\nStep4: Output the ribosome pausing score.', flush=True)
    pausing_score.output_cds_pausing()
    pausing_score.output_cds_codon_pausing()
    pausing_score.output_sum_codon_pausing()
    pausing_score.output_all_pausing()

    print('\nStep5: Draw the codon pausing score.', flush=True)
    pausing_score.draw_codon_total_pausing_heat()
    pausing_score.draw_codon_valid_pausing_heat()
    # pausing_score.draw_codon_pausing_plot('total')
    # pausing_score.draw_codon_pausing_plot('valid')

    if args.figure != 'none':
        print('\nStep6: Draw the RPFs pausing score.', flush=True)
        pausing_score.draw_rpf_pausing_plot()
    else:
        print('\nStep6: Do not draw figures.', flush=True)

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

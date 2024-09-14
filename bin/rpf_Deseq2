#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : DESeq2_pipe.py


from foo import ArgsParser
from foo import DESeq2_pipe


def main():
    ArgsParser.now_time()
    print('\nRun DESeq2 pipeline.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.deseq2_args_parser()
    deseq = DESeq2_pipe.DESeq2Pipe(args)

    print('\nStep2: Import the design file.', flush=True)
    deseq.import_design()

    print('\nStep3: Import the rpf file.', flush=True)
    deseq.import_reads()

    print('\nStep4: Prepare for DESeq2 analysis.', flush=True)
    deseq.perepare_deseq_pipe()

    print('\nStep5: Run DESeq2.', flush=True)
    deseq.run_deseq2()

    print('\nStep6: Output results.', flush=True)
    deseq.get_deseq_result()
    deseq.output_deseq()

    print('\nStep7: Draw volcano.', flush=True)
    deseq.prepare_volcano()
    deseq.draw_volcano()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()

if __name__ == '__main__':
    main()
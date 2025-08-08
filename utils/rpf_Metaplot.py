#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_metaplot.py


from utils.ribo import ArgsParser
from utils.ribo import Metaplot


def main():
    ArgsParser.now_time()
    print('\nDraw the metaplot.', flush=True)
    print('\nStep1: Checking the input Arguments.', flush=True)
    args = ArgsParser.metaplot_args_parser()
    meta = Metaplot.Metaplot(args)

    print('\nStep2: Import the RPFs file.', flush=True)
    meta.gene_anno()
    meta.read_rpf()
    meta.output_merge_meta()

    print('\nStep3: Draw the metaplot.', flush=True)
    if args.mode == "bar":
        meta.draw_bar_metaplot()
    elif args.mode == "line":
        meta.draw_line_metaplot()

    print('\nAll done.', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

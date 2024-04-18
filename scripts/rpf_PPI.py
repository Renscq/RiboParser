#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : rpf_PPI.py


from foo.ArgsParser import *
from foo.Pro_Pro_Interaction import *


def main():
    print('\nDraw the protein protein interaction network.\n', flush=True)
    print('\nStep1: Checking the input Arguments.\n', flush=True)
    args = ppi_args_parser()
    ppi = ProProInteraction(args)

    print('\nStep2: Set PPI version.\n', flush=True)
    ppi.set_ppi_version()

    print('\nStep3: Import the gene list.\n', flush=True)
    ppi.import_gene_list()

    print('\nStep4: Annotate the gene list.\n', flush=True)
    ppi.get_gene_id()

    print('\nStep5: Set the network parameters.\n', flush=True)
    ppi.set_ppi_version()
    ppi.set_network_parameter()

    print('\nStep6: Query for gene network.\n', flush=True)
    ppi.query_for_gene_network()

    print('\nStep7: Save the network graphic.\n', flush=True)
    ppi.save_network()

    print('\nAll done.\n', flush=True)
    now_time()


if __name__ == '__main__':
    main()


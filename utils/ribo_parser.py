#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : ribo_parser.py


from utils.ribo import ArgsParser
from utils.ribo.Ribo import *


def format_seq(ribo_attr):
    rpf_list = []

    for name, isoform in ribo_attr.mrna_dict.items():
        rpf_dict = OrderedDict()
        # print(name + '\n')

        # trim the cds frame to fit the 3nt periodicity
        if isoform.cds_length % 3 != 0:
            shift_cds = isoform.cds_length % 3
            isoform.cds_length -= shift_cds
            isoform.utr3_length += shift_cds

        shift_5 = isoform.utr5_length % 3
        shift_3 = isoform.utr3_length % 3

        utr5_length = isoform.utr5_length - shift_5
        utr3_length = isoform.utr3_length - shift_3
        trim_length = isoform.length - shift_5 - shift_3

        rpf_dict['name'] = [name] * int(trim_length / 3)
        rpf_dict['now_nt'] = list(range(1 + shift_5, isoform.length - shift_3, 3))
        rpf_dict['from_tis'] = list(range(-utr5_length // 3, (trim_length - utr5_length) // 3))
        rpf_dict['from_tts'] = list(range((utr3_length - trim_length) // 3 + 1, utr3_length // 3 + 1))
        rpf_dict['region'] = ['5utr'] * (utr5_length // 3) + ['cds'] * (isoform.cds_length // 3) + ['3utr'] * (utr3_length // 3)

        if shift_3 != 0:
            trim_seq = isoform.seq[shift_5: -shift_3]
            rpf_dict['aa_list'] = [str(trim_seq[i:i + 3]) for i in range(0, trim_length, 3)]

            trim_rpf = isoform.rpf[shift_5: -shift_3]
            rpf_dict[ribo_attr.output + '_f0'] = trim_rpf[0::3]
            rpf_dict[ribo_attr.output + '_f1'] = trim_rpf[1::3]
            rpf_dict[ribo_attr.output + '_f2'] = trim_rpf[2::3]
        elif shift_3 == 0:
            trim_seq = isoform.seq[shift_5::]
            rpf_dict['aa_list'] = [str(trim_seq[i:i + 3]) for i in range(0, trim_length, 3)]

            trim_rpf = isoform.rpf[shift_5::]
            rpf_dict[ribo_attr.output + '_f0'] = trim_rpf[0::3]
            rpf_dict[ribo_attr.output + '_f1'] = trim_rpf[1::3]
            rpf_dict[ribo_attr.output + '_f2'] = trim_rpf[2::3]
        try:
            rpf = pd.DataFrame(rpf_dict)
            rpf_list.append(rpf)
        except ValueError:
            print('Checking the input Arguments.\n'.format(gene_name=name), flush=True)

    total_rpf_df = pd.concat(rpf_list, ignore_index=True)
    # total_rpf_df["sum"] = total_rpf_df[['frame_1', 'frame_2', 'frame_3']].sum(axis=1)
    total_rpf_df.to_csv(ribo_attr.output + "_rpf.txt", sep='\t', index=False)


def main():
    ArgsParser.now_time()
    print('\nParse the RIBO-SEQ data.\n', flush=True)
    print('Step1: Checking the input Arguments.\n', flush=True)
    args = ArgsParser.ribo_args_parser()
    ribo_attr = Ribo(args)

    print('\nStep2: Import the p-site offset.\n', flush=True)
    ribo_attr.read_offset()

    print('\nStep3: Import the transcripts annotation.\n', flush=True)
    ribo_attr.read_transcript()

    print('\nStep4: Import the BAM file.\n', flush=True)
    ribo_attr.read_bam()

    print('\nStep5: Format the in-frame RPFs.\n', flush=True)
    format_seq(ribo_attr)

    print('\nAll done.\n', flush=True)
    ArgsParser.now_time()


if __name__ == '__main__':
    main()

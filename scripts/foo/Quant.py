#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Quant.py


import matplotlib.pyplot as plt
import seaborn as sns

import plotly.graph_objs as go
import plotly.express as px
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy

import numpy as np
import pandas as pd
from . import RPFs


class Quant(object):
    def __init__(self, args):
        self.rpf_file = args.rpf

        # self.anno = args.anno
        self.output_prefix = args.output

        # parameters for the trimmed CDS region
        self.rpf_num = 0
        self.site = 'P'
        self.tis = args.tis
        self.tts = args.tts
        self.frame = args.frame
        self.length = None

        # parameters for rpf file parser
        self.sample_name = None
        self.sample_num = None
        self.merged_rpf = None
        self.total_rpf = None
        self.col_index = None

        # parameters for different region
        self.utr5 = args.utr5
        self.utr3 = args.utr3

        self.utr5_rpf = None
        self.utr5_rpm = None
        self.utr5_tpm = None

        self.utr3_rpf = None
        self.utr3_rpm = None
        self.utr3_tpm = None

        self.cds_rpf = None
        self.cds_rpm = None
        self.cds_tpm = None

    def read_rpf(self):
        # raw_rpf, sample_name, sample_num, merged_rpf, total_rpf_num, gene_rpf_sum, high_gene, high_rpf
        rpf_results = RPFs.import_rpf(rpf_file=self.rpf_file,
                                      sites=self.site,
                                      frame=self.frame,
                                      sample_num=None,
                                      sample_name=None,
                                      tis=None,
                                      tts=None,
                                      gene=None,
                                      rpf_num=self.rpf_num)
        # self.raw_rpf = rpf_results[0]
        self.sample_name = rpf_results[1].copy()
        # self.sample_num = rpf_results[2]
        self.merged_rpf = rpf_results[3].copy()
        self.total_rpf = rpf_results[4]
        # self.gene_rpf_sum = rpf_results[5]
        # self.high_gene = rpf_results[6]
        # self.high_rpf = rpf_results[7].copy()
        del rpf_results

        self.sample_num = len(self.sample_name)
        self.col_index = ['name'] + self.sample_name
        self.length = self.merged_rpf.loc[:, ['name', 'region']].groupby("name").count() * 3

    def output_total_rpf(self):
        total_rpf_name = self.output_prefix + "_total.txt"
        self.total_rpf.index.name = 'sample'
        self.total_rpf.name = 'rpf_count'

        self.total_rpf.to_csv(total_rpf_name, sep='\t', index=True)

    def output_gene_rpf(self, rpf_count, rpf_rpm, rpf_rpkm, rpf_tpm, region):

        # output the rpf / rpm / tpm results
        rpf_count.astype(int)
        rpf_count.index.name = "name"
        rpf_count.to_csv(self.output_prefix + "_" + region + "_rpf_quant.txt", sep='\t', index=True)

        rpf_rpm = rpf_rpm.astype(float).round(2)
        rpf_rpm.index.name = "name"
        rpf_rpm.to_csv(self.output_prefix + "_" + region + "_rpm_quant.txt", sep='\t', index=True)

        rpf_rpkm = rpf_rpkm.astype(float).round(2)
        rpf_rpkm.index.name = "name"
        rpf_rpkm.to_csv(self.output_prefix + "_" + region + "_rpkm_quant.txt", sep='\t', index=True)

        rpf_tpm = rpf_tpm.astype(float).round(2)
        rpf_tpm.index.name = "name"
        rpf_tpm.to_csv(self.output_prefix + "_" + region + "_tpm_quant.txt", sep='\t', index=True)

    def quant_utr5_rpf(self, region='5utr'):

        print("Summary RPFs of 5'-UTR.", flush=True)
        # filter the gene region
        self.utr5_rpf = self.merged_rpf.loc[self.merged_rpf["region"] == region, self.col_index]
        self.utr5_rpf = self.utr5_rpf.groupby("name")[self.sample_name].sum()

        self.utr5_rpf.columns = [i + "_" + region + "_rpf" for i in self.sample_name]

        # calculate the gene rpm
        self.utr5_rpm = np.divide(self.utr5_rpf, np.asarray(self.total_rpf)) * 1e6
        self.utr5_rpm.columns = [i + "_" + region + "_rpm" for i in self.sample_name]

        # calculate the gene rpkm
        self.utr5_rpkm = (self.utr5_rpf * 1e9).div(self.length['region'], axis=0).div(np.asarray(self.total_rpf), axis=1)
        self.utr5_rpkm.columns = [i + "_" + region + "_rpkm" for i in self.sample_name]

        # calculate the gene tpm
        utr5_rpk = (self.utr5_rpf * 1000).div(self.length['region'], axis=0).dropna()
        utr5_rpk_sum = utr5_rpk.sum(axis=0)
        self.utr5_tpm = (utr5_rpk * 1000000).div(utr5_rpk_sum, axis=1)
        self.utr5_tpm.columns = [i + "_" + region + "_tpm" for i in self.sample_name]

        self.output_gene_rpf(self.utr5_rpf, self.utr5_rpm, self.utr5_rpkm, self.utr5_tpm, 'utr5')

    def quant_utr3_rpf(self, region='3utr'):

        print("Summary RPFs of 3'-UTR.", flush=True)
        # filter the gene region
        self.utr3_rpf = self.merged_rpf.loc[self.merged_rpf["region"] == region, self.col_index]
        self.utr3_rpf = self.utr3_rpf.groupby("name")[self.sample_name].sum()

        self.utr3_rpf.columns = [i + "_" + region + "_rpf" for i in self.sample_name]

        # calculate the gene rpm
        self.quant_rpm = np.divide(self.utr3_rpf, np.asarray(self.total_rpf)) * 1e6
        self.quant_rpm.columns = [i + "_" + region + "_rpm" for i in self.sample_name]

        # calculate the gene rpkm
        self.utr3_rpkm = (self.utr3_rpf * 1e9).div(self.length['region'], axis=0).div(np.asarray(self.total_rpf), axis=1)
        self.utr3_rpkm.columns = [i + "_" + region + "_rpkm" for i in self.sample_name]

        # calculate the gene tpm
        utr3_rpk = (self.utr3_rpf * 1000).div(self.length['region'], axis=0).dropna()
        utr3_rpk_sum = utr3_rpk.sum(axis=0)
        self.utr3_tpm = (utr3_rpk * 1000000).div(utr3_rpk_sum, axis=1)
        self.utr3_tpm.columns = [i + "_" + region + "_tpm" for i in self.sample_name]

        self.output_gene_rpf(self.utr3_rpf, self.utr3_rpm, self.utr3_rpkm, self.utr3_tpm, 'utr3')

    def quant_cds_rpf(self, region='cds'):

        print('Summary RPFs of CDS.', flush=True)
        # filter the gene region
        self.cds_rpf = self.merged_rpf[self.merged_rpf['from_tis'] >= self.tis]
        self.cds_rpf = self.cds_rpf[self.cds_rpf['from_tts'] < -self.tts]

        self.cds_rpf = self.cds_rpf.loc[self.cds_rpf["region"] == region, self.col_index]
        self.cds_rpf = self.cds_rpf.groupby("name")[self.sample_name].sum()

        self.cds_rpf.columns = [i + "_" + region + "_rpf" for i in self.sample_name]

        # calculate the gene rpm
        self.cds_rpm = np.divide(self.cds_rpf, np.asarray(self.total_rpf)) * 1e6
        self.cds_rpm.columns = [i + "_" + region + "_rpm" for i in self.sample_name]

        # calculate the gene rpkm
        # formula: RPKM = (10^9 * C) / (N * L)
        # C: number of reads mapped to a gene
        # N: total number of reads mapped to the reference
        # L: number of bases in a gene
        
        self.cds_rpkm = (self.cds_rpf * 1e9).div(self.length['region'], axis=0).div(np.asarray(self.total_rpf), axis=1)
        self.cds_rpkm.columns = [i + "_" + region + "_rpkm" for i in self.sample_name]

        # calculate the gene tpm
        cds_rpk = (self.cds_rpf * 1000).div(self.length['region'], axis=0).dropna()
        cds_rpk_sum = cds_rpk.sum(axis=0)
        self.cds_tpm = (cds_rpk * 1000000).div(cds_rpk_sum, axis=1)
        self.cds_tpm.columns = [i + "_" + region + "_tpm" for i in self.sample_name]

        self.output_gene_rpf(self.cds_rpf, self.cds_rpm, self.cds_rpkm, self.cds_tpm, 'cds')

    def quant_region(self):
        if self.utr3:
            self.quant_utr5_rpf(region='utr5')

        self.quant_cds_rpf(region='cds')

        if self.utr5:
            self.quant_utr3_rpf(region='utr3')

    def get_sub_plot_num(self):
        factor1 = int(self.sample_num ** 0.5)
        factor2 = int(self.sample_num ** 0.5)
        square = factor1 * factor2

        while not square >= self.sample_num:
            if factor2 > factor1:
                factor1 += 1
            else:
                factor2 += 1
            square = factor1 * factor2

        return factor1, factor2

    def draw_rpf_barplot(self):
        """
        draw the rpf region of each samples with plotly

        1. summary the rpf num of each region [5'-utr, cds, 3'-utr]
        2. draw the stacked bar plot
        """
        
        out_pdf = self.output_prefix + "_rpf_barplot.pdf"
        # out_png = self.output_prefix + "_rpf_barplot.png"

        rpf_sum = self.merged_rpf.groupby('region')[self.sample_name].apply(sum)
        rpf_sum_per = round(rpf_sum.div(rpf_sum.sum()) * 100, 2)

        rpf_sum_per_t = rpf_sum_per.reset_index().melt(id_vars="region",
                                                       var_name="samples",
                                                       value_name="proportion")

        fig = px.bar(rpf_sum_per_t,
                     x="samples",
                     y="proportion",
                     color="region",
                     orientation="v",
                     text_auto=True,
                     category_orders={"region": ["5utr", "cds", "3utr"]},
                     color_discrete_sequence=px.colors.qualitative.Pastel2)

        fig.update_layout(xaxis_title="proportion (%)",
                          yaxis_title="samples",
                          template="simple_white",
                          height=500, width=300 + len(self.sample_name) * 20)

        fig.write_image(out_pdf, engine="kaleido")
        # fig.write_image(out_png, engine="kaleido")

    def draw_rpf_cdfplot(self):
        """
        draw the cdf figure of each samples with plotly
        to check the global expression levels

        1. log the cds rpm
        2. draw the eCDF plot of each sample with plotly
        """
        out_pdf = self.output_prefix + "_rpf_cdfplot.pdf"
        # out_png = self.output_prefix + "_rpf_cdfplot.png"

        # get the log2 values
        cds_log = np.log2(self.cds_rpm + 1)
        cds_log = cds_log.reset_index().melt(id_vars="name", var_name="samples", value_name="rpm")

        # draw the proportion of RPFs in different region
        fig = px.ecdf(cds_log, x="rpm", color="samples")

        fig.update_layout(xaxis_title="eCDF",
                          yaxis_title="log2 RPM",
                          template="simple_white",
                          height=600, width=650)

        fig.write_image(out_pdf, engine="kaleido")
        # fig.write_image(out_png, engine="kaleido")

    def draw_rpf_pcaplot(self):
        """
        draw the pca figure of each samples with plotly

        1. log the cds rpm
        2. calculate the components of each dimension
        3. draw the PCA scatter plot of each sample
        """
        from sklearn.decomposition import PCA

        out_pdf = self.output_prefix + "_rpf_pcaplot.pdf"
        # out_png = self.output_prefix + "_rpf_pcaplot.png"

        # run the PCA with sklearn function
        out_txt = self.output_prefix + "_cds_rpf_pca.txt"

        pca = PCA(n_components=2)
        pca.fit(np.log2(self.cds_rpm.T + 1))
        variance_ratios = pca.explained_variance_ratio_
        # covariance = pca.get_covariance() * 100
        components = pca.fit_transform(np.log2(self.cds_rpm.T + 1))

        pca_df = pd.DataFrame(data=components, columns=['PC1', 'PC2'])
        pca_df.index = self.sample_name
        pca_df = pca_df.reset_index().rename(columns={'index': 'samples'})

        pca_df.to_csv(out_txt, sep='\t', index=True)

        pca_df_t = pca_df.reset_index().melt(id_vars="index",
                                             var_name="pc",
                                             value_name="value").rename(columns={'index': 'samples'})

        # draw the PCA figure
        fig = px.scatter(pca_df, x="PC1", y="PC2", color="samples")

        fig.update_xaxes(tickformat='.2e')
        fig.update_yaxes(tickformat='.2e')
        fig.update_layout(xaxis_title="PC1 (" + str("{:.2f}".format(variance_ratios[0])) + "%)",
                          yaxis_title="PC2 (" + str("{:.2f}".format(variance_ratios[1])) + "%)",
                          template="simple_white",
                          height=400, width=490)

        fig.write_image(out_pdf, engine="kaleido")

    def draw_rpf_heatmap1(self):
        """
        draw the gene rpm heatmap of each samples

        1. log the cds rpm
        2. scale the rpm to center 0
        3. draw the heatmap with di-cluster
        """

        out_pdf = self.output_prefix + "_rpf_heatmap.pdf"
        out_png = self.output_prefix + "_rpf_heatmap.png"

        cds_rpm = self.cds_rpm[(self.cds_rpm > 0).any(axis=1)]
        rpm_log = np.log2(cds_rpm + 1)

        plt.figure(figsize=(8, 12), dpi=300)

        sns.clustermap(rpm_log,
                       cmap="RdBu",
                       center=0,
                       dendrogram_ratio=(.1, .2),
                       cbar_pos=(.02, .02, .025, .17))

        plt.savefig(fname=out_pdf)
        plt.savefig(fname=out_png)

        plt.close()

    def draw_rpf_heatmap2(self):
        """
        draw the gene rpm heatmap of each samples with plotly

        1. log the cds rpm
        2. calculate the distance of rpm data matrix and clustered with scipy
        3. draw the di-cluster heatmap with plotly
        """

        out_pdf = self.output_prefix + "_rpf_heatmap.pdf"
        # out_png = self.output_prefix + "_rpf_heatmap.png"

        rpm_log = np.log2(self.cds_rpm[(self.cds_rpm > 0).any(axis=1)] + 1)

        cds_dist = pdist(rpm_log)
        clustering = hierarchy.linkage(cds_dist, method='ward', metric='euclidean')

        row_idx = hierarchy.leaves_list(clustering)
        rpm_clustered = rpm_log.iloc[row_idx, :]

        fig = go.Figure(
            data=go.Heatmap(x=rpm_clustered.columns.values,
                            y=rpm_clustered.index.values,
                            z=rpm_clustered,
                            # zmin=-3, zmax=3,
                            colorscale='RdBu',
                            reversescale=True)
        )

        fig.update_layout(xaxis=dict(side='bottom'),
                          yaxis=dict(side='left'),
                          height=600, width=450,
                          title=dict(text='Expression pattern', x=0.5, xanchor='center'),
                          xaxis_showgrid=False, yaxis_showgrid=False,
                          xaxis_tickangle=-90, yaxis_tickangle=0)

        fig.write_image(out_pdf, engine="kaleido")
        # fig.write_image(out_png, engine="kaleido")

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Project : riboParser
# @Script  : Pro_Pro_Interaction.py

import requests
import mygene
import pandas as pd


class ProProInteraction(object):

    def __init__(self, args):
        # get the input and output files
        self.input = args.input
        self.gene = []
        self.gene_anno = None
        self.number = args.number
        self.output = args.output
        self.figure = ""

        # set network parameter
        self.params = {}
        self.species = args.taxonomy
        self.scopes = args.scopes

        self.network_flavor = args.flavor
        self.required_score = args.score
        self.caller_identity = "www.awesome_app.org"
        self.show_query_node_labels = args.labels
        self.hide_disconnected_nodes = args.nodes
        self.block_structure_pics_in_bubbles = args.block
        self.network_type = args.type

        # set ppi version
        self.graphic_format = args.graphic
        self.string_api_url = "https://version-11-5.string-db.org/api"
        self.method = "network"
        self.request_url = ''
        self.response = ''

    def set_ppi_version(self):
        self.request_url = "/".join([self.string_api_url, self.graphic_format, self.method])

        # merge the output filename
        if self.graphic_format.upper() == 'SVG':
            self.figure = self.output + '.svg'
        else:
            self.figure = self.output + '.png'

    def import_gene_list(self):
        with open(self.input, 'r') as input_gene:
            now_num = 1
            for line in input_gene:
                if now_num > self.number:
                    break
                self.gene.append(line.strip())

    def get_gene_id(self):
        mg = mygene.MyGeneInfo()
        self.gene_anno = mg.querymany(self.gene,
                                     scopes=self.scopes,
                                     species=self.species,
                                     fields='name, taxid, symbol, entrezgene, refseq, uniprot, ensembl',
                                     as_dataframe=True)
        self.gene_anno.to_csv(self.output + '.txt', sep='\t')
        self.gene_anno = self.gene_anno['entrezgene'].astype('str').to_list()

    def set_network_parameter(self):
        self.params = {
            # your protein
            "identifiers": "\r".join(self.gene_anno),
            # species NCBI identifier
            "species": self.species,
            # show confidence links
            "network_flavor": self.network_flavor,
            # set the confidence score
            "required_score": self.required_score,
            # set the app name
            "caller_identity": self.caller_identity,
            # hides all protein names from the picture (0 or 1) (defailt:0)
            "show_query_node_labels": self.show_query_node_labels,
            # hides all proteins that are not connected to any other protein in your network (0 or 1) (default:0)
            "hide_disconnected_nodes": self.hide_disconnected_nodes,
            # disables structure pictures inside the bubble (0 or 1) (default:0)
            "block_structure_pics_in_bubbles": self.block_structure_pics_in_bubbles,
            # network type: functional (default), physical
            "network_type": self.network_type
        }

    def query_for_gene_network(self):
        self.response = requests.post(self.request_url, data=self.params)

    def save_network(self):
        print("Saving interaction network to %s" % self.figure, flush=True)

        with open(self.figure, 'wb') as network_graphic:
            network_graphic.write(self.response.content)

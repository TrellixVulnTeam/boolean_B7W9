import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import sys

from dataclasses import dataclass


import bone


@dataclass
class ALZanalysis:
    def __post_init__(self):
        self.gw_1 = self._get_json_weights(self._get_file("gene_weights_1.json"))

    def _get_file(self, file: str):
        file_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(file_dir, file)
        return file_path

    def _get_json_weights(self, json_file: str):
        with open(json_file) as file_in:
            gene_weights = json.load(file_in)
        gene_weights = {int(k): v for k, v in gene_weights.items()}
        return gene_weights

    def dong2013(self):
        gse40060 = bone.GEO(accessionID="GSE40060")
        survival = gse40060.survival()
        expr = gse40060.expr(rename_genes=True, probeID="ENSG")
        expr = expr.fillna(0)
        dong = bone.BoNE(expr, survival)

        name = "c source_name_ch1"
        groups = ["endogenous", "overexpressed"]
        dong.init(name, self.gw_1, groups)
        return dong

    def dong2013_2(self):
        gse40060 = bone.GEO(accessionID="GSE40060")


def violin(bone_obj):
    plt.figure(figsize=(10, 5), dpi=100)
    bone_obj.violin()
    plt.savefig("dong_violin.png")


if __name__ == "__main__":
    alz = ALZanalysis()
    alz.dong2013_2()
    # d = alz.dong2013()
    # violin(d)

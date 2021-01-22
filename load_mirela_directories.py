import yaml
import os 
from singlecellexperiment import SingleCellExperiment
import sys
import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy
import glob
from singlecellexperiment import SingleCellExperiment
from scipy import sparse
from scipy import io
import collections
import subprocess
import scanpy as sc


sces = []
mats = []
metadata = collections.defaultdict(list)
genes = []
metadata["clone"] = []

sces = glob.glob("*/*/*.rdata")
clones = glob.glob("*/*/*.csv")

def load_dict(samples):
    meta = dict()
    for s in samples:
        res = s.split("/")[-1].split(".")[0]
        meta[res] = s
    return meta

sces = load_dict(sces)
clones = load_dict(clones)

samples = set(sces.keys()).intersection(set(clones.keys()))
yamls = glob.glob("sample_*.yaml")
for yamlfile in yamls:
    with open(yamlfile) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        for sample, smetadata in data.items():
            if sample not in samples:
                print("Skipping {}.".format(sample))
                continue
            print("Loading {}".format(sample))
            obj = SingleCellExperiment.fromRData(sces[sample])
            mat = obj.assays["logcounts"]
            mats.append(mat)
            df = pandas.read_csv(clones[sample])
            metadata["sample"] += [sample for _ in range(mat.shape[1])]
            metadata["timepoint"] += [smetadata["timepoint"] for _ in range(mat.shape[1])]
            metadata["series"] += [smetadata["series"] for _ in range(mat.shape[1])]
            metadata["condition"] += [smetadata["treat"] for _ in range(mat.shape[1])]
            metadata["material"] += [smetadata["material"] for _ in range(mat.shape[1])]
            metadata["tenx"] += [smetadata["tenx"] for _ in range(mat.shape[1])]
            metadata["barcode"] += obj.colData["Barcode"]
            
            clone_dict = dict(zip(df["Barcode"],df["clone"]))
            for barcode in obj.colData["Barcode"]:
                if barcode in clone_dict:
                    metadata["clone"].append(clone_dict[barcode])
                else:
                    metadata["clone"].append("None")
            if len(genes) == 0:
                genes = obj.rownames
                symbols = obj.rowData["ensembl_gene_id"]



df = pandas.DataFrame.from_dict(metadata)
df.to_csv("metadata.csv")
print("metadata done.")

genemap = dict(zip(symbols,genes))

matrix = sparse.hstack(mats).tocsr()

print(matrix.shape)
for column, values in metadata.items():
    print(column, len(values))
print(len(genes))

mtx = "result/matrix.mtx"
barcodesgz = "result/barcodes.tsv"
featuresgz = "result/features.tsv"

io.mmwrite(mtx, matrix)
subprocess.call(['gzip', mtx])

output = open(barcodesgz,"w")
for barcode in metadata["barcode"]:
    output.write(barcode+"\n")
output.close()
subprocess.call(['gzip', barcodesgz])

output = open(featuresgz,"w")
for gene in symbols:
    output.write(gene + "\t" + genemap[gene] + "\tGene Expression\n")
output.close()
subprocess.call(['gzip', featuresgz])
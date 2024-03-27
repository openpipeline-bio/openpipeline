import json
import os
import mudata as mu
from pathlib import Path
from typing import Dict
import warnings
import torch
import numpy as np
from torch.utils.data import Dataset, DataLoader
from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

#TODO: Check preprocessing requirements
# data_is_raw=False, filter_gene_by_counts=False, filter_cell_by_counts=False, subset_hvg=False
# always binning?
# alwasy cross-check genes?
#TODO: check 'hyperparameter_path' (resulting in parameters) and its values
# 'embsize', 'd_hid', 'nlayers', 'nhead', 'dropout' get overridden by model_config
# need contents for default values of seed, amp
#TODO: DSBN True or False? 
# Currently True, but all batch labels are removed so is essentially False
#TODO: parametrization of tokenize module
# 'max_seq_len' and 'include_zero_gene
#TODO: refactor code in case gene name layer is an index

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_preprocessed.h5mu",
    "modality": "rna",
    "input_gene_ids": 'resources_test/scgpt/test_resources/Kim2020_Lung_gene_ids.pt',
    "input_values": 'resources_test/scgpt/test_resources/Kim2020_Lung_values.pt',
    "input_padding_mask": 'resources_test/scgpt/test_resources/Kim2020_Lung_padding_mask.pt',
    "model": "resources_test/scgpt/source/best_model.pt",
    "model_config": "resources_test/scgpt/source/args.json",
    "model_vocab": "resources_test/scgpt/source/vocab.json",
    "output": "Kim2020_Lung_cell_type_annotation.h5ad",
    "gene_name_layer": "gene_name",
    "batch_id_layer": "batch_id",
    "predicted_cell_type_id": "predicted_cell_type",
    "pad_token": "<pad>",
    "mask_ratio": 0,
    "mask_value": -1,
    "pad_value": -2,
    "n_cls": 8,
    "n_input_bins": 51,
    "batch_size": 64, # Check with config
    "output_compression": None
}
## VIASH END

# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}

warnings.filterwarnings('ignore')

# set_seed(config.seed)
set_seed(0)

logger.info("Reading in data")

# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

#TODO: Optionally set device to use gpu, to be implemented when gpu-based machine types are supported
# if par["device"] == "cuda" and not torch.cuda.is_available():
#     print("WARNING: CUDA is not available. Using CPU instead.")
#     device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device(par["device"])
logger.info("Setting device to use cpu")
device = torch.device("cpu")

# Set padding specs
special_tokens = [par["pad_token"], "<cls>", "<eoc>"]

# Set batch labels to 1
adata.obs["str_batch"] = "1"
batch_id_labels = adata.obs['str_batch'].astype("category").cat.codes.values
adata.obs["batch_id"] = batch_id_labels
batch_ids = adata.obs["batch_id"].tolist()
num_batch_types = len(set(batch_ids))

# Load vocabulary
vocab_file = par["model_vocab"]
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)
genes = adata.var[par["gene_name_layer"]].tolist()
vocab.set_default_index(vocab["<pad>"])
gene_ids = np.array(vocab(genes), dtype=int)
ntokens = len(vocab) 

# load model configs
model_config_file = par["model_config"]
with open(model_config_file, "r") as f:
    model_configs = json.load(f)
logger.info(
    f"Model args in {model_config_file} will override the config."
)
embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]
n_layers_cls = model_configs["n_layers_cls"]
dropout = model_configs["dropout"]

# Instantiate model
model = TransformerModel(
    ntokens,
    d_model=embsize,
    nhead=nhead,
    d_hid=d_hid,
    nlayers=nlayers,
    nlayers_cls=3,
    n_cls=par["n_cls"],
    vocab=vocab,
    dropout=dropout,
    pad_token=par["pad_token"],
    pad_value=par["pad_value"],
    do_mvc=False, # changed from config to hard-coded
    do_dab=False, # hard-coded
    use_batch_labels=False, 
    num_batch_labels=num_batch_types,
    domain_spec_batchnorm=False, # Check in combination with 
    input_emb_style="continuous", # hard-coded
    n_input_bins=par["n_input_bins"],
    cell_emb_style="cls", # hard-coded
    mvc_decoder_style="inner product", # needed? MVC set to false
    ecs_threshold=0, # needed? ECS set to false in forward method
    explicit_zero_prob=False, # True for integration embedding - check out why
    use_fast_transformer=False,  #TODO: parametrize when GPU is available
    fast_transformer_backend="flash", #TODO: parametrize when GPU is available
    pre_norm=False, #TODO: parametrize when GPU is available
)

model_file = par["model"]
try:
    model.load_state_dict(torch.load(model_file, map_location=device))
    logger.info(f"Loading all model params from {model_file}")
except:
    # only load params that are in the model and match the size
    model_dict = model.state_dict()
    pretrained_dict = torch.load(model_file, map_location=device)
    pretrained_dict = {
        k: v
        for k, v in pretrained_dict.items()
        if k in model_dict and v.shape == model_dict[k].shape
    }
    for k, v in pretrained_dict.items():
        logger.info(f"Loading params {k} with shape {v.shape}")
    model_dict.update(pretrained_dict)
    model.load_state_dict(model_dict)

model.to(device)

batch_ids = adata.obs["batch_id"].tolist()
batch_ids = np.array(batch_ids)

input_gene_ids = torch.load(par["input_gene_ids"])
input_values = torch.load(par["input_values"])

test_data_pt = {
    "gene_ids": input_gene_ids,
    "values": input_values,
    "batch_labels": torch.from_numpy(batch_ids).long(),
}

test_loader = DataLoader(
    dataset=SeqDataset(test_data_pt),
    batch_size=par["batch_size"],
    shuffle=False,
    drop_last=False,
    # num_workers=min(os.cpu_count(), par["batch_size"] // 2),
    # num_workers=min(len(os.sched_getaffinity(0)), par["batch_size"] // 2),
    pin_memory=True,
)

model.eval()

predictions = []
with torch.no_grad():
    for batch_data in test_loader:
        input_gene_ids = batch_data["gene_ids"].to(device)
        input_values = batch_data["values"].to(device)
        batch_labels = batch_data["batch_labels"].to(device)

        src_key_padding_mask = input_gene_ids.eq(vocab[par["pad_token"]])
        with torch.cuda.amp.autocast(enabled=False): # parametrize
            output_dict = model(
                input_gene_ids,
                input_values,
                src_key_padding_mask=src_key_padding_mask,
                batch_labels=None,
                CLS=True,  # evaluation does not need CLS or CCE
                CCE=False,
                MVC=False,
                ECS=False,
                do_sample=False, # changed from configurable (config) to hard-coded
                #generative_training = False,
            )
            output_values = output_dict["cls_output"]
        preds = output_values.argmax(1).cpu().numpy()
        predictions.append(preds)
        
predictions = np.concatenate(predictions, axis=0)

# Write output
logger.info("Writing output data")
adata.obs[par["predicted_cell_type_id"]] = predictions
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])

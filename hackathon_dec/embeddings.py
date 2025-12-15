# --------------------------------------------------
# 5️⃣ ESM Embedding (Sequence → Features)
# --------------------------------------------------
from transformers import AutoTokenizer, AutoModel
import esm
import os
import pandas as pd
import torch
import gemmi
from Bio import SeqIO
def get_esm_embedding(sequence: str):
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()

    batch_converter = alphabet.get_batch_converter()
    data = [("protein", sequence)]
    _, _, tokens = batch_converter(data)

    with torch.no_grad():
        out = model(tokens, repr_layers=[33])
        token_emb = out["representations"][33]

    embedding = token_emb[0, 1:-1].mean(dim=0)
    print(embedding)
    return embedding


# --------------------------------------------------
# 6️⃣ ProtBERT Embedding (Sequence → Features)
# --------------------------------------------------
def get_protbert_embedding(sequence: str):
    tokenizer = AutoTokenizer.from_pretrained(
        "Rostlab/prot_bert_bfd",
        do_lower_case=False
    )
    model = AutoModel.from_pretrained("Rostlab/prot_bert_bfd")
    model.eval()

    seq = " ".join(list(sequence))
    inputs = tokenizer(seq, return_tensors="pt")

    with torch.no_grad():
        outputs = model(**inputs)
        hidden = outputs.last_hidden_state

    embedding = hidden[0, 1:-1].mean(dim=0)
    print(embedding)
    return embedding

get_protbert_embedding("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEW")
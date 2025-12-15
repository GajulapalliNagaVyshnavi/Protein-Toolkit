import os
import pandas as pd
import torch
import gemmi
from Bio import SeqIO
from transformers import AutoTokenizer, AutoModel
import esm


# --------------------------------------------------
# 1️⃣ FASTA → CSV
# --------------------------------------------------
def fasta_to_csv(fasta_path: str, csv_path: str | None = None):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(fasta_path)

    if csv_path is None:
        base, _ = os.path.splitext(fasta_path)
        csv_path = base + ".csv"

    records = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        records.append({
            "seq_id": record.id,
            "sequence": str(record.seq)
        })

    pd.DataFrame(records).to_csv(csv_path, index=False)
    print(f"FASTA → CSV saved: {csv_path}")


# --------------------------------------------------
# 2️⃣ CSV → FASTA
# --------------------------------------------------
def csv_to_fasta(csv_path: str, target_column=None,fasta_path: str | None = None):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(csv_path)
    if target_column==None:
        target_column="sequence"

    if fasta_path is None:
        base, _ = os.path.splitext(csv_path)
        fasta_path = base + ".fa"

    df = pd.read_csv(csv_path)

    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            f.write(f">{row['seq_id']}\n{row['target_column']}\n")

    print(f"CSV → FASTA saved: {fasta_path}")


# --------------------------------------------------
# 3️⃣ PDB → CIF
# --------------------------------------------------
def pdb_to_cif(pdb_path: str, cif_path: str | None = None):
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(pdb_path)

    if cif_path is None:
        base, _ = os.path.splitext(pdb_path)
        cif_path = base + ".cif"

    structure = gemmi.read_structure(pdb_path)
    structure.make_mmcif_document().write_file(cif_path)
    print(f"PDB → CIF saved: {cif_path}")


# --------------------------------------------------
# 4️⃣ CIF → PDB
# --------------------------------------------------
def cif_to_pdb(cif_path: str, pdb_path: str | None = None):
    if not os.path.exists(cif_path):
        raise FileNotFoundError(cif_path)

    if pdb_path is None:
        base, _ = os.path.splitext(cif_path)
        pdb_path = base + ".pdb"

    structure = gemmi.read_structure(cif_path)
    structure.write_pdb(pdb_path)
    print(f"CIF → PDB saved: {pdb_path}")

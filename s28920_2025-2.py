#!/usr/bin/env python3
"""
NCBI GenBank Extended Retriever
"""

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

def fetch_records(email, api_key, taxid, min_len, max_len, max_records=200):
    Entrez.email = email
    Entrez.api_key = api_key

    search_term = f"txid{taxid}[Organism]"
    handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y", retmax=max_records)
    search = Entrez.read(handle)
    ids = search["IdList"]
    webenv = search["WebEnv"]
    query_key = search["QueryKey"]

    handle = Entrez.efetch(
        db="nucleotide", rettype="gb", retmode="text",
        webenv=webenv, query_key=query_key, retmax=max_records
    )

    records = SeqIO.parse(handle, "genbank")
    filtered = []

    for record in records:
        length = len(record.seq)
        if min_len <= length <= max_len:
            filtered.append({
                "Accession": record.id,
                "Length": length,
                "Description": record.description
            })

    return filtered

def save_csv(data, filename="report.csv"):
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"[+] CSV report saved as {filename}")

def plot_lengths(data, filename="length_chart.png"):
    sorted_data = sorted(data, key=lambda x: x["Length"], reverse=True)
    accs = [x["Accession"] for x in sorted_data]
    lens = [x["Length"] for x in sorted_data]

    plt.figure(figsize=(10, 5))
    plt.plot(accs, lens, marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("GenBank Accession")
    plt.ylabel("Sequence Length")
    plt.title("Sequence Lengths by Accession")
    plt.tight_layout()
    plt.savefig(filename)
    print(f"[+] Length chart saved as {filename}")

def main():
    email = input("Enter your email: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxonomic ID (taxid): ")
    min_len = int(input("Minimum sequence length: "))
    max_len = int(input("Maximum sequence length: "))

    print("[*] Fetching records...")
    records = fetch_records(email, api_key, taxid, min_len, max_len)

    if not records:
        print("[-] No records found in specified length range.")
        return

    save_csv(records)
    plot_lengths(records)
    print("[✓] All done.")

if __name__ == "__main__":
    main()

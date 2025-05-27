#!/usr/bin/env python3

from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt


def get_user_input():
    email = input("Email: ")
    api_key = input("API Key: ")
    tax_id = input("TaxID: ")

    try:
        min_len = int(input("Min Length: "))
        max_len = int(input("Max Length: "))
    except ValueError:
        exit("Invalid input.")

    return email, api_key, tax_id, min_len, max_len


def search_sequences(tax_id, min_len, max_len):
    query = f"txid{tax_id}[Organism] AND {min_len}:{max_len}[SLEN]"
    handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y", retmax=0)
    results = Entrez.read(handle)
    handle.close()

    if int(results["Count"]) == 0:
        exit("No records found.")

    return results["WebEnv"], results["QueryKey"]


def fetch_sequences(webenv, query_key, min_len, max_len):
    handle = Entrez.efetch(
        db="nucleotide", rettype="gb", retmode="text",
        retmax=100, webenv=webenv, query_key=query_key
    )

    records = []
    for record in SeqIO.parse(handle, "genbank"):
        seq_len = len(record.seq)
        if min_len <= seq_len <= max_len:
            records.append({
                "Accession": record.id,
                "Length": seq_len,
                "Description": record.description
            })

    handle.close()
    return records


def save_results(data, tax_id):
    df = pd.DataFrame(data).sort_values("Length", ascending=False)
    csv_name = f"taxid_{tax_id}_report.csv"
    df.to_csv(csv_name, index=False)
    print(f"Saved data to {csv_name}")

    # Plotting
    plt.figure(figsize=(10, 5))
    plt.plot(df["Accession"], df["Length"], marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Length")
    plt.title(f"Seq Lengths for TaxID {tax_id}")
    plt.tight_layout()
    chart_name = f"taxid_{tax_id}_chart.png"
    plt.savefig(chart_name)
    print(f"Saved chart to {chart_name}")
    return len(df)


def main():
    email, api_key, tax_id, min_len, max_len = get_user_input()
    Entrez.email = email
    Entrez.api_key = api_key

    webenv, query_key = search_sequences(tax_id, min_len, max_len)
    records = fetch_sequences(webenv, query_key, min_len, max_len)
    count = save_results(records, tax_id)
    print(f"Saved {count} records to CSV and PNG.")


if __name__ == "__main__":
    main()

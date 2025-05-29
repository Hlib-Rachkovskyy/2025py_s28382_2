#!/usr/bin/env python3
from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email, Entrez.api_key, Entrez.tool = email, api_key, 'BioScriptEx10'

    def search_taxid(self, taxid):
        print(f"Searching for records with taxID: {taxid}")
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            print(f"Found {count} records for {organism_name} (TaxID: {taxid})")
            if count == 0:
                return None


            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=100):
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            batch_size = min(max_records, 500)
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            records = list(SeqIO.parse(handle, "gb"))
            return records

        except Exception as e:
            print(f"Error fetching records: {e}")
            return []

    def filter_by_length(self, records, min_len, max_len):
        return [rec for rec in records if min_len <= len(rec.seq) <= max_len]

    def generate_csv(self, records, filename):
        data = [(record.id, len(record.seq), record.description) for record in records]
        pd.DataFrame(data, columns=["Accession number", "Length", "Description"]).to_csv(filename, index=False)
        print(f"Saved CSV report to {filename}")

    def generate_plot(self, records, filename):
        sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=True)
        args = [rec.id for rec in sorted_records]
        lengths = [len(rec.seq) for rec in sorted_records]

        plt.figure(figsize=(20, 5))
        plt.plot(args, lengths, marker='o')
        plt.xticks(rotation=90)
        plt.title("GenBank Record Lengths")
        plt.xlabel("Accession number")
        plt.ylabel("Sequence length")
        plt.tight_layout()
        plt.savefig(filename)
        print(f"Saved plot to {filename}")
        plt.close()


def main():
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    retriever = NCBIRetriever(email, api_key)

    taxid = input("Enter taxonomic ID (taxid) of the organism: ")

    min_len = int(input("Enter minimum sequence length: "))
    max_len = int(input("Enter maximum sequence length: "))

    count = retriever.search_taxid(taxid)
    if not count:
        print("No records found. Exiting.")
        return

    print("\nFetching records...")
    all_records = retriever.fetch_records(start=0, max_records=200)

    print(f"Fetched {len(all_records)} records. Filtering by length...")
    filtered = retriever.filter_by_length(all_records, min_len, max_len)

    if not filtered:
        print("No records matched the length criteria.")
        return

    csv_file = f"taxid_{taxid}_report.csv"
    retriever.generate_csv(filtered, csv_file)

    plot_file = f"taxid_{taxid}_plot.png"
    retriever.generate_plot(filtered, plot_file)


if __name__ == "__main__":
    main()

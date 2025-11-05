#!/usr/bin/env python3
import argparse
import time
from pathlib import Path
from Bio import Entrez
import os
import requests
import gzip
import shutil
import random
import pandas as pd
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf
import json
import mygene
import certifi
import ssl

os.environ.setdefault('SSL_CERT_FILE', certifi.where())
# -----------------------
# Utilidades
# -----------------------
def read_genes_from_file(filename):
    """
    Lee una lista de genes desde un fichero.
    Cada línea puede contener uno o varios genes separados por comas.
    """
    with open(filename, "r") as file:
        content = file.read().strip()
        # Divide por comas y limpia espacios en blanco
        genes = [g.strip() for g in content.split(",") if g.strip()]
    return genes

# -----------------------
# Funciones
# -----------------------
def bio_python(gene_list, email: str = "testing@uma.es"):
    """
    input_file: Ruta al fichero con la lista de genes
    email: correo requerido por NCBI Entrez
    """
    output_file = os.path.join("results", "gene_annotations.txt")

    Entrez.email = email

    with open(output_file, "w") as out:
        for gene in gene_list:
            try:
                print(f"Searching for gene: {gene}...")
                handle = Entrez.esearch(db="gene", term=f"{gene}[Gene] AND Homo sapiens[Organism]")
                record = Entrez.read(handle)
                handle.close()

                if record['IdList']:
                    gene_id = record['IdList'][0]

                    print(f"Fetching annotations for {gene} (Gene ID: {gene_id})...")
                    handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
                    gene_record = handle.read()
                    handle.close()

                    print(f"Annotations for {gene}:\n")
                    print(gene_record)
                    print("\n" + "="*50 + "\n")

                    out.write(f"Annotations for {gene}:\n")
                    out.write(gene_record)
                    out.write("\n" + "="*50 + "\n")

                    print(f"Annotations for {gene} saved.\n")
                else:
                    print(f"No record found for gene: {gene}")
                    out.write(f"No record found for gene: {gene}\n")

                # Pausa entre consultas para respetar los servidores de Entrez
                time.sleep(1)

            except Exception as e:
                print(f"An error occurred while fetching annotations for {gene}: {e}")
                out.write(f"Error for {gene}: {e}\n")

    print(f"\nAll results saved to '{output_file}'.")

def goaTools(gene_list, fdr_value=0.05):
    """
    Ejecuta un análisis de enriquecimiento GO usando GOATOOLS.
    gene_list : lista de símbolos génicos (ej. ['COX4I2','ND1','ATP6'])
    fdr_value : umbral de FDR para filtrar resultados/significancia (float)
    """

    out_csv = os.path.join("results", "go_enrichment_results.csv")
    out_txt = os.path.join("results", "go_enrichment.txt")

    # --- Descarga / comprobación de archivos GO y GAF ---
    obo_url = "http://current.geneontology.org/ontology/go-basic.obo"
    obo_file_path = "go-basic.obo"
    if not os.path.exists(obo_file_path):
        try:
            print("Downloading go-basic.obo ...")
            response = requests.get(obo_url)
            response.raise_for_status()
            with open(obo_file_path, "wb") as fh:
                fh.write(response.content)
            print(f"Downloaded GO OBO to {obo_file_path}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download go-basic.obo: {e}")
            return
    else:
        print(f"GO OBO file '{obo_file_path}' already exists. Skipping download.")

    gaf_url = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
    compressed_gaf = "goa_human.gaf.gz"
    extracted_gaf = "goa_human.gaf"
    if not os.path.exists(extracted_gaf):
        try:
            print("Downloading goa_human.gaf.gz ...")
            response = requests.get(gaf_url, stream=True)
            response.raise_for_status()
            with open(compressed_gaf, "wb") as fh:
                for chunk in response.iter_content(chunk_size=8192):
                    fh.write(chunk)
            print(f"Downloaded GAF gzip to {compressed_gaf}")

            print("Extracting GAF file ...")
            with gzip.open(compressed_gaf, "rb") as f_in:
                with open(extracted_gaf, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"Extracted GAF to {extracted_gaf}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download or extract GAF file: {e}")
            return
    else:
        print(f"GAF file '{extracted_gaf}' already exists. Skipping download/extraction.")

    # --- Map gene symbols to UniProt (Swiss-Prot) using MyGene ---
    print("Mapping gene symbols to UniProt IDs using MyGene.info ...")
    mg = mygene.MyGeneInfo()
    # querymany puede lanzar excepciones por problemas de red; lo hacemos robusto con retry básico
    try:
        query_result = mg.querymany(gene_list, scopes='symbol', fields='uniprot', species='human')
    except Exception as e:
        print(f"Error querying MyGene.info: {e}")
        return

    study_gene_uniprot_ids = []
    # Manejar distintos formatos de 'uniprot' retornado por MyGene
    for entry in query_result:
        # entry puede contener .get('uniprot') como dict o lista; intentamos extraer Swiss-Prot
        uni = entry.get('uniprot')
        uniprot_id = None
        if isinstance(uni, dict):
            # puede ser 'Swiss-Prot' con string o lista
            swiss = uni.get('Swiss-Prot')
            if isinstance(swiss, list) and swiss:
                uniprot_id = swiss[0]
            elif isinstance(swiss, str):
                uniprot_id = swiss
        elif isinstance(uni, list) and uni:
            # lista de dicts o strings; intentamos tomar primer valor
            first = uni[0]
            if isinstance(first, dict):
                uniprot_id = first.get('Swiss-Prot') or first.get('TrEMBL') or None
            elif isinstance(first, str):
                uniprot_id = first

        if uniprot_id:
            study_gene_uniprot_ids.append(uniprot_id)

    if not study_gene_uniprot_ids:
        print("No UniProt IDs could be retrieved for the provided gene symbols.")
        return

    print("Converted UniProt IDs for the study genes:")
    for gs, up in zip(gene_list, study_gene_uniprot_ids):
        print(f"Gene Symbol: {gs}, UniProt ID: {up}")

    # --- Construir background a partir del GAF ---
    all_genes = set()
    try:
        with open(extracted_gaf, "r") as gf:
            for line in gf:
                if line.startswith("!"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) > 1:
                    gid = cols[1]
                    if gid:
                        all_genes.add(gid)
        print(f"Extracted {len(all_genes)} gene identifiers from GAF file.")
    except FileNotFoundError:
        print(f"GAF file '{extracted_gaf}' not found.")
        return
    except Exception as e:
        print(f"Error reading GAF file: {e}")
        return

    all_genes_list = list(all_genes)
    # Si hay menos de 1000 genes en total, usamos todos; si no, muestreamos 1000
    random.seed(42)
    if len(all_genes_list) <= 1000:
        background_gene_list = all_genes_list.copy()
    else:
        background_gene_list = random.sample(all_genes_list, 1000)

    # Asegurarnos de que los genes de estudio están en el background
    for gid in study_gene_uniprot_ids:
        if gid not in background_gene_list:
            background_gene_list.append(gid)

    print(f"Background list prepared with {len(background_gene_list)} genes (including study genes).")

    # --- Cargar ontología y asociaciones ---
    print("Loading GO ontology (obo) ...")
    try:
        go = obo_parser.GODag(obo_file_path)
    except Exception as e:
        print(f"Failed to parse OBO file: {e}")
        return

    print("Reading gene2go associations from GAF ...")
    try:
        gene2go = read_gaf(extracted_gaf)
    except Exception as e:
        print(f"Failed to read GAF associations: {e}")
        return

    # --- Ejecutar GO enrichment ---
    print("Running GO enrichment analysis with GOATOOLS ...")
    go_enrich = GOEnrichmentStudy(
        background_gene_list,
        gene2go,
        go,
        propagate_counts=False,
        alpha=0.05,
        methods=['fdr_bh']
    )

    # Ejecutar el estudio (puede tardar dependiendo del tamaño)
    try:
        enriched_results = go_enrich.run_study(study_gene_uniprot_ids)
    except Exception as e:
        print(f"GO enrichment study failed: {e}")
        return

    # --- Filtrar y guardar resultados ---
    rows = []
    txt_lines = []
    txt_lines.append(f"GO enrichment results (FDR threshold = {fdr_value})\n")
    txt_lines.append("=" * 80 + "\n")

    for res in enriched_results:
        # Algunos registros pueden no tener p_fdr_bh; manejar con getattr y None por defecto
        fdr = getattr(res, "p_fdr_bh", None)
        if fdr is None:
            # intentar p_fdr_bh puede no existir si no se calcularon correcciones
            continue
        if fdr < fdr_value:
            study_items = getattr(res, "study_items", [])
            genes_in_go = ", ".join(study_items)
            row = {
                "GO": getattr(res, "GO", ""),
                "Name": getattr(res, "name", ""),
                "FDR": fdr,
                "P_uncorrected": getattr(res, "p_uncorrected", None),
                "Study_count": getattr(res, "study_count", None),
                "Pop_count": getattr(res, "pop_count", None),
                "Study_items": genes_in_go
            }
            rows.append(row)
            txt_lines.append(f"GO ID: {row['GO']}, Description: {row['Name']}, FDR: {row['FDR']}, Genes: {genes_in_go}\n")

    # Guardar CSV
    if rows:
        df = pd.DataFrame(rows)
        try:
            df.to_csv(out_csv, index=False)
            print(f"Enrichment results saved to CSV: {out_csv}")
        except Exception as e:
            print(f"Failed to write CSV results: {e}")
    else:
        print("No enriched GO terms passed the FDR threshold. No CSV file created.")

    # Guardar TXT legible
    try:
        with open(out_txt, "w") as fh:
            fh.writelines(txt_lines)
        print(f"Enrichment summary saved to text file: {out_txt}")
    except Exception as e:
        print(f"Failed to write text summary: {e}")

    # Imprimir los resultados filtrados por pantalla (en inglés)
    if rows:
        print("\nSignificant GO terms (filtered by FDR):")
        for r in rows:
            print(f"GO ID: {r['GO']}, Description: {r['Name']}, FDR: {r['FDR']}, Genes: {r['Study_items']}")
    else:
        print("No significant GO terms to display.")

    # Final
    print("\nGO enrichment analysis completed.")
    # breve pausa para evitar saturar APIs si hay llamadas posteriores en el pipeline
    time.sleep(0.5)

def stringDB(proportioned_genes, fdr_value=0.05):
    """
    Realiza un análisis de enriquecimiento funcional utilizando la API de STRINGdb
    para una lista de genes dada, y guarda los resultados en results/stringdb_results.csv.
    """

    # URL de la API de STRINGdb y detalles del método
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment"

    # Construir el URL de la solicitud
    request_url = "/".join([string_api_url, output_format, method])

    # Genes para usar en el análisis de enriquecimiento funcional
    genes = proportioned_genes

    # Definir los parámetros para la solicitud
    params = {
        "identifiers": "%0d".join(genes),  # Lista de proteínas formateada con %0d (separador de nueva línea)
        "species": 9606,                   # Homo sapiens
        "caller_identity": "test_HAB"      # Identificador de prueba o correo
    }

    print(f"[INFO] Querying STRINGdb for {len(genes)} genes...\n")

    # Realizar la solicitud POST a la API de STRINGdb
    response = requests.post(request_url, data=params)

    if response.status_code != 200:
        print(f"[ERROR] STRINGdb request failed. Status code: {response.status_code}")
        return

    # Analizar la respuesta JSON
    data = json.loads(response.text)

    if not data:
        print("[WARNING] No data received from STRINGdb.")
        return

    # Filtrar resultados y mostrar por pantalla los términos significativos
    print(f"\n[RESULTS] Functional Enrichment Analysis (GO Biological Process, FDR < {fdr_value}):\n")
    results = []
    for row in data:
        term = row.get("term")
        preferred_names = ",".join(row.get("preferredNames", []))
        fdr = float(row.get("fdr", 1))
        description = row.get("description", "")
        category = row.get("category", "")

        if category == "Process" and fdr < fdr_value:
            print(f"{term:<12} {preferred_names:<40} {fdr:<10.4g} {category:<10} {description}")
            results.append({
                "term": term,
                "genes": preferred_names,
                "FDR": fdr,
                "category": category,
                "description": description
            })

    # Guardar resultados en CSV
    os.makedirs("results", exist_ok=True)
    out_csv = os.path.join("results", "stringdb_results.csv")

    if results:
        df = pd.DataFrame(results)
        df.to_csv(out_csv, index=False)
        print(f"\n[INFO] Results saved to: {out_csv}")
    else:
        pd.DataFrame(columns=["term", "genes", "FDR", "category", "description"]).to_csv(out_csv, index=False)
        print(f"\n[INFO] No significant results found (FDR < {fdr_value}). Empty CSV created at: {out_csv}")

# -----------------------
# CLI
# -----------------------
def cli():
    p = argparse.ArgumentParser(description="Fetch gene annotations from NCBI using Biopython Entrez")
    p.add_argument("input_file", type=Path, help="Path to input file with gene symbols (one or many per line).")
    p.add_argument("--email", type=str, default="your_email@example.com", help="Email to set for Entrez (required by NCBI).")
    args = p.parse_args()
    
    if not args.input_file.exists():
        print(f"Input file '{args.input_file}' not found.")
        return
    
    # --- Preparar carpeta results ---
    os.makedirs("results", exist_ok=True)

    genes = read_genes_from_file(args.input_file)
    print(f"Genes read: {genes}")
    if genes:
        bio_python(genes, email=args.email)
        goaTools(genes, fdr_value=0.05)
        stringDB(genes, fdr_value=0.05)

if __name__ == "__main__":
    cli()
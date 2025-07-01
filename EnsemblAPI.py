import ensembl_rest
import pandas as pd

'''
Workflow Summary

# 1. fetch orthology for a gene (get_orthologs)
# 2. for each id fouund, fetch gene info (get_gene_structure)
# 4. 
# 5. create a dataframe
'''

def get_gene_info(symbol, species):

    '''Args: 
    Species - string specifying species
    Symbol - string specifying gene name (e.g. ACTB for beta-actin)
    '''

    gene_info = ensembl_rest.symbol_lookup(species = species, symbol = symbol)
    return gene_info

def get_orthologs(gene_id, target_species):
    
    try:
        gene_orthologies = ensembl_rest.homology_symbol(
                        species='homo_sapiens',
                        symbol=gene_id,
                        params={'target_species': target_species,
                        'sequence': 'none',
                        'type': 'orthologues',
                        'format':'condensed'
                        }
        )
        print(gene_orthologies)
        return gene_orthologies

    # Handle the case when there's no gene tree
    except ensembl_rest.HTTPError as err:
        error_code = err.response.status_code
        error_message = err.response.json()['error']
        if (error_code == 400) \
        and ('Lookup found nothing' in error_message):
            # Skip the gene with no data
            pass
        else:
            # The exception was caused by another problem
            # Raise the exception again
            raise

def get_orthologs_info(homology_response):

    orthologs = []

    for entry in homology_response['data']:
        for hom in entry['homologies']:
            species = hom['species']
            id = hom['id']
            orthologs.append((id, species))

    return orthologs

def get_gene_structure(gene_id, species):
    """
    Use lookup_id with expand and condensed format to fetch transcripts and exons.
    """
    return ensembl_rest.lookup(
        gene_id,
        params={
            'species': species,
            'expand': True,
            'format': 'full'
        }
    )

def get_iso_exo_count(gene_data):

    transcripts = gene_data.get('Transcript', [])
    iso_count = len(transcripts)

    canonical_id = gene_data.get('canonical_transcript', '').split(".")[0]

    canonical = next(
        (tr for tr in transcripts if tr.get('id') == canonical_id or tr.get('is_canonical') == 1),
        None
    )

    if canonical:
        exo_count = len(list(canonical.get('Exon', [])))
        return iso_count, exo_count, canonical["id"]
    else:
        iso_count, None, None



def main():

    gene_id_list = ['ARMH1', 'ACTB', 'DSCAM']
    target_species = ['Homo sapiens', 'Pongo abelii', 'Mus musculus', 'Pan troglodytes', 'Gallus gallus', 'Canis lupus familiaris', 'Danio rerio']
    table_rows = []

    with open("ortholog_gene_info.txt", "w") as ortholog_file, open("IsoExoCount.txt", 'w') as isoexofile: 

        for gene_id in gene_id_list:

            # Add human separately
            human_gene_info = get_gene_info(gene_id, 'homo_sapiens')
            human_gene_id = human_gene_info['id']
            human_gene_data = get_gene_structure(human_gene_id, 'homo_sapiens')

            ortholog_file.write(f"homo_sapiens - {human_gene_id}\n\n")
            ortholog_file.write(str(human_gene_data) + "\n\n")

            num_isoforms, num_exons, canonical_id = get_iso_exo_count(human_gene_data)
            isoexofile.write(f"homo_sapiens - {human_gene_id}\n\n")
            isoexofile.write(f"N of isoforms for gene {gene_id}: {num_isoforms}, N of exons in canonical transcript {canonical_id}: {num_exons} \n\n")

            table_rows.append({
                "gene_name": gene_id,
                "species": "homo_sapiens",
                "ortholog_gene_id": human_gene_id,
                "num_transcripts": num_isoforms,
                "canonical_transcript_id": canonical_id,
                "num_exons_canonical": num_exons
            })

            #All others
            orthology = get_orthologs(gene_id, target_species=target_species)
            orthologs = get_orthologs_info(orthology)

            for ortholog in orthologs:
                gene_data = get_gene_structure(ortholog[0], ortholog[1])
                ortholog_file.write(f"{ortholog[1]} - {ortholog[0]}\n\n")  # species + gene ID
                ortholog_file.write(str(gene_data) + "\n\n")

                num_isoforms, num_exons, canonical_id = get_iso_exo_count(gene_data)
                isoexofile.write(f"{ortholog[1]} - {ortholog[0]}\n\n")  # species + gene ID
                isoexofile.write(f"N of isoforms for gene {gene_id}: {num_isoforms}, N of exons in canonical transcript {canonical_id}: {num_exons} \n\n")

                table_rows.append({
                        "gene_name": gene_id,
                        "species": ortholog[1],
                        "ortholog_gene_id": ortholog[0],
                        "num_transcripts": num_isoforms,
                        "canonical_transcript_id": canonical_id,
                        "num_exons_canonical": num_exons
                    })
    

    df = pd.DataFrame(table_rows)
    print(df)
    df.to_csv("ortholog_exon_summary.csv", index=False)

    

if __name__ == '__main__':
    main()
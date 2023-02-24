from biorun.libs import placlib as plac
import requests, json, csv, time
import pickle
import pandas as pd

# API: https://maayanlab.cloud/Enrichr/help#api


@plac.opt('counts', "input counts", abbrev='c')
@plac.opt('organism', "input counts", abbrev='d')
@plac.opt('colname', "gene id column name", abbrev='n')
@plac.opt('pval_cutoff', "pvalue cutoff", abbrev='t', type=float)
@plac.opt('pval_column', "pvalue column name", abbrev='p')
@plac.opt('output', "pvalue column name", abbrev='o')
def run(counts="edger.csv", organism='mmusculus', colname='gene', pval_cutoff=0.05, pval_column='FDR',
        output='enrichr.csv'):
    """
    Runs the enrichr tool on a csv file where one column contains gene names and some column contains pvalues.

    Filters the p values by a threshold, then submits the gene names to g:GOSt.

    Valid organisms: https://biit.cs.ut.ee/gprofiler/page/organism-list
    """

    ct = pd.read_csv(counts)

    ct = ct[ct[pval_column] < pval_cutoff]

    genes = ct[colname].tolist()

    def strip_dot(x):
        return x.split('.')[0]

    # Get rid of version numbers if these exists
    genes = list(map(strip_dot, genes))

    print(f"# Running Enrichr")
    # print(f"# https://biit.cs.ut.ee/gprofiler/gost")
    print(f"# Counts: {counts}")
    print(f"# Organism: {organism}")

    print(f"# Name column: {colname}")
    print(f"# Pval column: {pval_column} < {pval_cutoff}")
    print(f"# Gene count: {len(genes)}")
    print(f"# Genes: {','.join(genes[:5])},[...]")

    def submit():

        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
        genes_str = '\n'.join(genes)
        description = 'Example gene list'
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }


        print(f"# Submitting to Enrichr")

        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        data = json.loads(response.text)

        print(f"# User list id: {data.get('userListId')}")

        return data

    def export(data, library='KEGG_2015'):
        """
        Exports the results from enrichr.
        """
        # Get the user list id.
        userListId = data['userListId']

        # Construct the URL.
        export_url = f'https://maayanlab.cloud/Enrichr/export?userListId={userListId}&backgroundType={library}'

        # Submit the request.
        response = requests.get(export_url, stream=True)

        # Read the response.
        text = response.content.decode('utf-8')

        # Turn the response into comma separated file.
        lines = text.splitlines()
        lines = map(lambda x: x.strip(), lines)
        lines = map(lambda x: x.split('\t'), lines)

        def reformat(x):
            if "/" in x[1]:
                x[1] = x[1].replace("/", " of ")
            return x
        lines = map(reformat, lines)
        lines = list(lines)


        # Make text a UTF-8 string.
        fp = open(output, 'wt')
        writer = csv.writer(fp, quoting=csv.QUOTE_MINIMAL)



        writer.writerows(lines)
        fp.close()

        print(f"# Entries: {len(lines)} ")
        print(f"# Output: {output}")

    data = submit()
    export(data)

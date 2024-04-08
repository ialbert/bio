import requests
from xml.etree.ElementTree import fromstring as xml_from_string
from biorun.libs import placlib as plac
from biorun.libs import xmltodict
import json

# http://www.biomart.org/martservice_9.html

HOST = 'https://www.biomart.org'
PATH = '/biomart/martservice'

url = 'http://www.ensembl.org/biomart/martservice'

query2 = '''
<!DOCTYPE Query>
<Query client="true" processor="TSV" limit="-1" header="1">
    <Dataset name="drenio_gene_ensembl" config="gene_ensembl_config">
        <Filter name="chromosome_name" value="1" filter_list=""/>
        <Attribute name="ensembl_transcript_id"/>
        <Attribute name="ensembl_gene_id"/>
    </Dataset>
</Query>
'''


def showmarts():



    #params=dict(type='datasets', mart=self._name)

    params = dict(type='registry')


    resp = requests.get(url, params=params)

    text = resp.content

    out = xmltodict.parse(text)

    out = json.dumps(out, indent=4)
    print(out)
    pass


@plac.opt('x', "clade to remove", abbrev='R')
def run(x):

    #showmarts()

    query = '''
    <!DOCTYPE Query>
    <Query client="true" processor="TSV" limit="-1" header="1">
        <Dataset name="drerio_gene_ensembl" config="gene_ensembl_config">
            <Filter name="chromosome_name" value="1" filter_list=""/>
            <Attribute name="ensembl_transcript_id"/>
            <Attribute name="ensembl_gene_id"/>
        </Dataset>
    </Query>
    '''

    params = dict(query=query)

    resp = requests.get(url, params=params, stream=True)

    for chunk in resp.iter_content(decode_unicode=True):
        print (chunk, end='')

    #text = resp.content



if __name__ == '__main__':
    plac.call(run)

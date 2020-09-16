Demo accession numbers:

Viral, SARS-COV-2: NC_045512 (1 second)

    time bio fetch NC_045512 -v > foo.gb

Viral, Ebola: AF086833 (1 second)

    time bio fetch AF086833 -v > foo.gb

Bacteria, Escherichia coli: NC_002695 (2 minutes, 10MB)

    time bio fetch NC_002695 -v > foo.gb

Insect, fruit fly (chromosome 2L): NT_033779 (2 minutes, 44MB)

    time bio fetch NT_033779  -v > foo.gb

Mammal, human genome (chromosome 1): NC_000001 ( hours, 300MB)

    time bio fetch NC_000001  -v > foo.gb
    

E-utilites:

web site:


* https://www.ncbi.nlm.nih.gov/books/NBK25497/?report=reader

code example:

    URL=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?
    
    PARAMS='db=nuccore&id=AF086833&rettype=gbwithparts&retmode=text'
    
    #PARAMS1='db=nuccore&id=AF086833&rettype=gb&retmode=text'
    
    echo time curl ${URL}${PARAMs}
    
    
data_mode

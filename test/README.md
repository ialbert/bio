Demo accession numbers:

Viral, SARS-COV-2: NC_045512 (1 second)

    time bio fetch NC_045512

Viral, Ebola: AF086833 (1 second)

    time bio fetch AF086833 

Bacteria, Escherichia coli: NC_002695 (2 minutes, 10MB)

    time bio fetch NC_002695 -v > foo.gb

Insect, fruit fly (chromosome 2L): NT_033779 (2 minutes, 44MB)

    time bio fetch NT_033779 

Mammal, human genome (chromosome 1): NC_000001 (2 minutes, 300MB)

    time bio fetch NC_000001  -v > foo.gb
    

https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=gbwithparts&id=459277393&extrafeat=null&withparts=on

E-utilites:

web site:


* https://www.ncbi.nlm.nih.gov/books/NBK25497/?report=reader

code example:

    URL=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?
    
    PARAMS='db=nuccore&id=AF086833&rettype=gbwithparts&retmode=text'
    
    #PARAMS1='db=nuccore&id=AF086833&rettype=gb&retmode=text'
    
    echo time curl ${URL}${PARAMs}
    
    
data_mode

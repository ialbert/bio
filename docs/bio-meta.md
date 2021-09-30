# bio meta: obtain viral metadata

TODO: Work in progress.

Display metadata on a taxonomical id (with a header):

    bio meta 11138 -H | head

prints:

    accession,species,host,date,location,isolate,species_name
    NC_048217.1,11138,,,,,Murine hepatitis virus
    NC_001846.1,11138,,,,,Murine hepatitis virus
    AC_000192.1,11144,,,,,Murine hepatitis virus strain JHM
    JQ173883.1,1163669,,,,,Murine hepatitis virus strain S/3239-17
    AY700211.1,11138,,,,,Murine hepatitis virus
    AF208067.1,11138,,,,,Murine hepatitis virus
    AF208066.1,11138,,,,,Murine hepatitis virus
    AF201929.1,76344,,,,,Murine hepatitis virus strain 2
    AF029248.1,11138,,,,,Murine hepatitis virus


The column order is:

* accession
* species
* host,date
* location
* isolate
* species_name


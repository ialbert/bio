# `bio align`: perform alignments {#bio-align}

**Note:** alignments in `bio` are primarily designed for exploratory use, for aligning relatively short (up to ~30Kb long sequences), visually investigating the alignments, interacting with the sequences before and after alignment.

Use software that relies on heuristics when investigating large datasets. Specialized software will operate (many) orders of magnitude faster. Depending on your needs you may want to use: `blast`, `blat`, `mummer`, `minimap2`, `lastz`, `lastal`, `exonerate`, `vsearch`, `diamond`.


## Quick start

Input may be given from command line:

    bio align GATTACA GATCA

prints:

    # DNA: SeqA (7) vs SeqB (5) score=24.0
    # Alignment: length=7 ident=4/7(57.1%) mis=1 del=2 ins=0 gap=2
    # Parameters: match=5 mismatch=4 gap-open=11 gap-extend=1

    GATTACA
    |||.|   7
    GATCA--

output may be formatted in tabular fashion or as variants

    bio align GATTACA GATCA --variant

prints:

    pos type   len target  query
    4    mis    1    T       C
    6    del    2    CA      --

the default alignment is semi-global (global alignment with no end gap penalies). To perform a global alignment write:

    bio align GATTACA GATCA --global

the output is:

    # DNA: SeqA (7) vs SeqB (5) score=13.0
    # Alignment: length=7 ident=5/7(71.4%) mis=0 del=2 ins=0 gap=2
    # Parameters: match=5 mismatch=4 gap-open=11 gap-extend=1

    GATTACA
    |||  || 7
    GAT--CA

## Align realistic data

Fetches two files

    bio fetch NC_045512 MN996532 > genomes.gb

Aligning two genomes:

    cat genomes.gb | bio fasta | bio align | head

the command aligns two 30KB sequences and takes about 15 seconds on my system, it will print:

    # DNA: NC_045512.2 (29,903) vs MN996532.2 (29,855) score=148070.0
    # Alignment: length=29903 ident=28720/29903(96.0%) mis=1135 del=48 ins=0 gap=48
    # Parameters: match=5 mismatch=4 gap-open=11 gap-extend=1

    ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAA
    ||||||||||||||||||.|||||||||||||||||.||||.||||||||||||||||||||||||||||||||||||||| 81
    ATTAAAGGTTTATACCTTTCCAGGTAACAAACCAACGAACTCTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAA

The `bio align` method takes the first file that it sees as target and aligns all other sequences to it as queries.

## Proteins alignments

Align the DNA corresponding to protein `S`

    cat genomes.gb | bio fasta --gene S | bio align | head

prints:

    # DNA: YP_009724390.1 (3,822) vs QHR63300.2 (3,810) score=18767.0
    # Alignment: length=3822 ident=3549/3822(92.9%) mis=261 del=12 ins=0 gap=12
    # Parameters: match=5 mismatch=4 gap-open=11 gap-extend=1

    ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTACCCCCTGCA
    ||||||||||||||||||||||||||||||||.||||||||||||||||||||.|||||.||||||||.|||||.|||||| 81
    ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTTTCTAGTCAGTGTGTTAATCTAACAACTAGAACTCAGTTACCTCCTGCA

or as translated sequences that prints:

    cat genomes.gb | bio fasta --gene S --translate | bio align | head

that prints:

    # PEP: YP_009724390.1 (1,274) vs QHR63300.2 (1,270) score=6542.0
    # Alignment: length=1274 ident=1241/1274(97.4%) mis=29 del=4 ins=0 gap=4
    # Parameters: matrix=BLOSUM62 gap-open=11 gap-extend=1

    MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDN
    |||||||||||||||||||||||||||||||.|||||||||||||||||.|||||||||||||||||||||||||.||||| 81
    MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSSTRGVYYPDKVFRSSVLHLTQDLFLPFFSNVTWFHAIHVSGTNGIKRFDN

## Alignment with tabular output

    cat genomes.gb | bio fasta --gene S --translate | bio align --table

prints (check the default output to identify which number corresonds to what):

    YP_009724390.1	QHR63300.2	1274	97.4	29	4	0

## Alignment with variant output

    cat genomes.gb | bio fasta --gene S --translate | bio align --variant

prints the variation (some lines not show):

    pos   type    len    target    query
    32     mis     1       F         S
    50     mis     1       S         L
    346    mis     1       R         T
    372    mis     1       A         T
    403    mis     1       R         T
    439    mis     3       NNL       KHI
    604    mis     1       T         A
    681    del     4       PRRA      ----
    1125   mis     1       N         S
    1228   mis     1       V         I
    ...

## Different scoring matrices

    cat genomes.gb | bio fasta --gene S --translate | bio align --matrix PAM30 | head

prints:

    # PEP: YP_009724390.1 (1,274) vs QHR63300.2 (1,270) score=9339.0
    # Alignment: length=1274 ident=1241/1274(97.4%) mis=29 del=4 ins=0 gap=4
    # Parameters: matrix=PAM30 gap-open=11 gap-extend=1

    MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDN
    |||||||||||||||||||||||||||||||.|||||||||||||||||.|||||||||||||||||||||||||.||||| 81
    MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSSTRGVYYPDKVFRSSVLHLTQDLFLPFFSNVTWFHAIHVSGTNGIKRFDN

## Scoring matrices

The scoring matrix may be a builtin name or a file with a scoring matrix format. See the scoring with:

    bio align --matrix PAM30

prints:

    #
    # This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
    #
    # PAM 30 substitution matrix, scale = ln(2)/2 = 0.346574
    #
    # Expected score = -5.06, Entropy = 2.57 bits
    #
    # Lowest score = -17, Highest score = 13
    #
        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    A   6  -7  -4  -3  -6  -4  -2  -2  -7  -5  -6  -7  -5  -8  -2   0  -1 -13  -8  -2  -3  -3  -3 -17
    R  -7   8  -6 -10  -8  -2  -9  -9  -2  -5  -8   0  -4  -9  -4  -3  -6  -2 -10  -8  -7  -4  -6 -17
    N  -4  -6   8   2 -11  -3  -2  -3   0  -5  -7  -1  -9  -9  -6   0  -2  -8  -4  -8   6  -3  -3 -17
    D  -3 -10   2   8 -14  -2   2  -3  -4  -7 -12  -4 -11 -15  -8  -4  -5 -15 -11  -8   6   1  -5 -17
    ...



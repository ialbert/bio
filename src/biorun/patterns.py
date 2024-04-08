"""
Stores the patterns for various data resources
"""
import re

# Match for SRR numbers: SRR5260547
SRR_PATT = re.compile(r'(ERR|SRR|DRR|SRP|ERP)\d+')

# function that matches SRR numbers
is_srr = lambda text: bool(SRR_PATT.fullmatch(text))

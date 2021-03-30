import os
import sys


def parse(fname):
    stream = open(fname)
    stream = map(lambda x: x.rstrip(), stream)
    stream = filter(lambda x: x.strip(), stream)
    text = '\n'.join(stream)

    return text


cmd = sys.argv[1]

os.system(f"bash {cmd} > tempfile")

inp = parse(cmd)
out = parse("tempfile")

start = f"""
```bash
{inp}
```
"""

end = f"""
the command above prints:

```
{out}
```
"""

if out:
    print(start + end)
else:
    print(start)

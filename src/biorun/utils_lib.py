import json
import csv
import sys

def as_json(reader, indent=4, lines=False, stream=sys.stdout):
    """
    Returns a JSON string or JSON lines.
    """
    
    if lines:
        for item in reader:
            stream.write(json.dumps(item))
            stream.write("\n")
    else:
        data = list(reader)
        stream.write(json.dumps(data, indent=indent))


def as_csv(reader, stream=sys.stdout, sep=','):
    """
    Writes a CSV stream from a reader.
    """
    writer = csv.DictWriter(stream, fieldnames=reader.fieldnames, delimiter=sep)
    writer.writeheader()
    writer.writerows(reader)

if __name__ == '__main__':
    
    as_json(123)
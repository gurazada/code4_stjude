"""
A simple parser for the GFF3 format.

Based on Uli Koehler's version.
Updated by CFT to:
    handle blank lines
    cleaned up class names
    improved string output
Updated by GG to:
    urllib2
    preserve attributes order
    handle empty attributes in string output
    mutable record attributes

Test with transcripts.gff3 from
http://www.broadinstitute.org/annotation/gebo/help/gff3.html.

Format specification source:
http://www.sequenceontology.org/gff3.shtml

Version 1.0
"""
from collections import OrderedDict
from recordclass import recordclass
import gzip
import urllib2

#Initialized GeneInfo record class.
#Note: since namedtuple is immutable, using recordclass instead. Attributes are now mutable.
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

class GFF3Record(recordclass("GFF3Record", gffInfoFields)):
    def __str__(self):
        """Provide a more useful string output"""
        return "\t".join(('.' if self.seqid is None else self.seqid,
                          '.' if self.source is None else self.source,
                          '.' if self.type is None else self.type,
                          '.' if self.start is None else str(self.start),
                          '.' if self.end is None else str(self.end),
                          '.' if self.score is None else str(self.score),
                          '.' if self.strand is None else self.strand,
                          '.' if self.phase is None else str(self.phase),
                          '.' if not self.attributes else ';'.join(["{}={}".format(k, v) for k, v in self.attributes.items()])))


def parse_GFF_attributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""#
    if attributeString == ".": return {}
    ret = OrderedDict()
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[urllib2.unquote(key)] = urllib2.unquote(value)
    return ret

def parse_GFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith("#") or not line: continue
            parts = line.split("\t")
            #If this fails, the file format is not standard-compatible
            #print('len parts = {}, len fileds = {}'.format(len(parts), len(gffInfoFields)))
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else urllib2.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib2.unquote(parts[1]),
                "type": None if parts[2] == "." else urllib2.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib2.unquote(parts[6]),
                "phase": None if parts[7] == "." else int(parts[7]),
                "attributes": parse_GFF_attributes(parts[8])
            }
            yield GFF3Record(**normalizedInfo)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="The GFF3 input file (.gz allowed)")
    parser.add_argument("--print-records", action="store_true", help="Print all GeneInfo objects, not only")
    args = parser.parse_args()
    #Execute the parser
    recordCount = 0
    parser = parse_GFF3(args.file)
    for record in parser:
        #Print record if specified by the user
        if args.print_records:
            print(record)
            #Access attributes like this: my_strand = record.strand
        recordCount += 1

    print("Total records: {}".format(recordCount))

if __name__ == "__main__":
    main()
# -- coding: utf-8 --

# by Guna on 10/2015

# import some modules
import os
import sys
import argparse
import textwrap
import re
import json
import urllib2
import gzip
from operator import itemgetter
from collections import OrderedDict
# other Pioneer developed modules and constants
import phi.Logger
from .gff3parser import parse_GFF3
from .DN_Constants import DN_ANALYSIS_NAME, \
    GAIA_ASSEMBLY_API, SPECIESCODE_TAXON_API, \
    SPECIESCODE_SCINAME_API, GENE_COUNTER_START, \
    GENE_COUNTER_STEP, DP_PREFIX, SOURCE_COLUMN, \
    GFF_HEADER


# this read the json data into a dict and convert unicode to ascii
# convert array data
# copied from internet.
def byteify(input):
    if isinstance(input, dict):
        return {byteify(key):byteify(value) for key,value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input


class DisplayName(object):
    """class for creating display names for genes and mRNA/tRNA features in gff3"""

    def __init__(
            self,
            inputFile='',
            outputFile='',
            gsapRunId='',
            seqIdRegex='',
            assemblyId=None,
            taxonId=None,
            scientificName=None,
            logger=None
        ):
        
        self.inputFile = inputFile # input gff3 file
        self.outputFile = outputFile # output gff3 file with display names created
        self.gsapRunId = gsapRunId
        self.seqIdRegex = seqIdRegex # regex in the format 'ZmChr(?P<num>\d+)v2' for ZmChr01v2
        # parse seqid func matching the regex => group('num') gives the chr/contig number
        self.parse_seqid = re.compile(r'{}'.format(seqIdRegex)).match
        # One of the 3: assemblyId, taxonId, sciname is required to get species code
        # look up Genome Repo based on assemblyId, to get taxonId
        self.assemblyId = assemblyId
        # Else, either internal taxon ID or scientific name is required to get species code
        self.taxonId = taxonId
        self.scientificName = scientificName
        self.logger = logger
        
        
        self.speciesCode = ''
        # data structures to create
        # array of all GFF3 records
        self.records = []
        # dict of all seqid-->chr# in gff3
        self.seqidDict = {}
        # temp bin to hold non-matching sequence IDs
        self.nomatch_seqid_bin = []
        # array of all gene features: {seqid, start, id, name, children=[{id, name}], records=[]}
        # children only upto transcript-level mRNA or tRNA
        self.geneFeatureArr = []
        # final id-->displayname map
        self.id2nameMap = OrderedDict()
        
        self.outputDir = os.path.dirname(self.outputFile)
        self.id2nameOutputFile = self.outputFile + '_id2nameMap'
        self.seqidOutputFile = self.outputFile + '_seqidMap'
        
    
    def check(self):
        """check inputs"""
        # check on input file - required
        if not self.inputFile or not os.path.exists(self.inputFile):
            self.logger.error("Input GFF3 file not set or does not exist: %s" % self.inputFile)
            raise Exception("Input GFF3 file not set or does not exist: %s" % self.inputFile)
        self.inputFile = os.path.abspath(self.inputFile)
        
        # check on outputDir
        if not os.path.exists(self.outputDir):
            self.logger.error("Output dir does not exist: %s" % self.outputDir)
            raise Exception("Output dir does not exist: %s" % self.outputDir)

        if not os.access(self.outputDir, os.W_OK) or not os.access(self.outputDir, os.X_OK):
            self.logger.error("Output dir %s is either not writable or executable. \
I need both to be able to write to it" % self.outputDir)
            raise Exception("Output dir %s is either not writable or executable. \
I need both to be able to write to it" % self.outputDir)
        
        # check on ids required for species code
        if not self.assemblyId and not self.taxonId and not self.scientificName:
            self.logger.error("Need either: assembly_id, taxon_id or scientific_name to get species_code for display name")
            raise Exception("Need either: assembly_id, taxon_id or scientific_name to get species_code for display name")
        
        # get species code
        if self.assemblyId is not None:
            self.speciesCode = self.speciescode_by_assemblyid(self.assemblyId)
            self.logger.info("Retrieved species code: {0} using assembly id: {1}".format(self.speciesCode, self.assemblyId))
        elif self.taxonId is not None:
            self.speciesCode = self.speciescode_by_taxonid(self.taxonId)
            self.logger.info("Retrieved species code: {0} using taxon id: {1}".format(self.speciesCode, self.taxonId))
        elif self.scientificName is not None:
            self.speciesCode = self.speciescode_by_sciname(self.scientificName)
            self.logger.info("Retrieved species code: {0} using scientific name: {1}".format(self.speciesCode, self.scientificName))


    def parseGff3FixDupIDs(self):
        """parse input GFF3 file and uniquify any duplicate IDs for gene, mRNA/tRNA features
           Important Assumption: gene structure is in order and grouped, parents coming first"""
           
        self.check()
        self.logger.info("Parsing input gff3 file and checking for duplicate IDs in gene, mRNA features: %s..." % self.inputFile)
        records = parse_GFF3(self.inputFile)
        # serializing entire GFF3 file as I/O is more expensive than RAM memory
        self.records = list(records)
        uniqueIdSet = {}  # gene, mRNA ID set - dictionary {ID: True} 
        renameIdMap = {}  # dictionary of IDs that have been renamed - {Old ID: new ID}
        dupCount = 0
        for record in self.records:
            # checking duplicate IDs
            if record.type in ['gene', 'mRNA', 'tRNA']:
                # only these features are main parents
                if 'ID' in record.attributes:
                    if record.attributes['ID'] in uniqueIdSet:
                        # found duplicate ID
                        dupCount += 1
                        newId = record.attributes['ID'] + '.' + str(dupCount) + 'd'
                        renameIdMap[record.attributes['ID']] = newId
                        self.logger.warn("Duplicate ID: {0} found. Assigning new ID: {1}".format(record.attributes['ID'], newId))
                        record.attributes['ID'] = newId
                        uniqueIdSet[record.attributes['ID']] = True
                    else:
                        # ID is unique - add to uniqueIdSet
                        uniqueIdSet[record.attributes['ID']] = True
                else:
                    self.logger.error("ID attribute missing in gff3 record: {}".format(record))
                    raise Exception("ID attribute missing in gff3 record: {}".format(record))
            else:
                # other child features follow and will be made unique
                # assumed format: parentId.exon1, parentId.cds1, parentId.utr5p1, parentId.utr3p1 etc
                if 'ID' in record.attributes:
                    parts = record.attributes['ID'].split('.')
                    parentId = '.'.join(parts[:-1])
                    if parentId in renameIdMap:
                        # parentId was renamed - update child ID accordingly
                        record.attributes['ID'] = '.'.join([renameIdMap[parentId], parts[-1]])
                else:
                    self.logger.error("ID attribute missing in gff3 record: {}".format(record))
                    raise Exception("ID attribute missing in gff3 record: {}".format(record))
                
            # Updating Parents if present - checks every record    
            if 'Parent' in record.attributes:
                # Parent attribute could have multiple CSVs
                new_parents = []
                for parent in record.attributes['Parent'].split(','):
                    if parent in renameIdMap:
                        new_parents.append(renameIdMap[parent])
                    else:
                        new_parents.append(parent)
                record.attributes['Parent'] = ','.join(new_parents)
                        
        self.logger.info("Parsed a total of {0} records. Found and fixed {1} duplicate IDs in gene/mRNA features.".format(len(self.records), dupCount))    
        
        
        
    def buildGeneSet(self):
        """create gene feature array"""
        self.parseGff3FixDupIDs()
        self.logger.info("Building gene features dict...")
        for record in self.records:
            if record.type == 'gene':
                # create seqidDict - parsing only gene features will be enough
                # parse and store the seqid<-->num if not already in seqidDict
                seqid = record.seqid
                if (seqid not in self.seqidDict) and (seqid not in self.nomatch_seqid_bin):
                    # new seqid
                    seqid_match = self.parse_seqid(seqid)
                    if seqid_match and 'num' in seqid_match.groupdict():
                        # match - store in seqidDict
                        self.seqidDict[seqid] = int(seqid_match.group('num'))
                    else:
                        # no match - we will get to them later
                        self.logger.warn("given regex '{0}' does not help extract number from seqid '{1}'. Please match (?P<num>\d+) group with number part if possible. Assigning automatic numbering, check seqidMap file: {2}".format(self.seqIdRegex, seqid, self.seqidOutputFile))
                        self.nomatch_seqid_bin.append(seqid)
                    
                # create geneDict
                # expecting validated GFF files - every gene feature has unique ID attribute
                geneDict = {'seqid': record.seqid,
                            'start': record.start,
                            'ID': record.attributes['ID'],
                            'Name': '', # display name
                            'children': [], # mRNA/tRNA features
                            'records': [] # all gff3 records comprising a gene
                            }
                # adding gene feature record
                geneDict['records'].append(record)
                self.geneFeatureArr.append(geneDict)
            elif record.type in ['mRNA', 'tRNA']:
                if not 'Parent' in record.attributes:
                    self.logger.error("found gff3 feature with missing Parent attribute. Need parents for child features: {0}\n".format(record))
                    raise Exception("found gff3 feature with missing Parent attribute. Need parents for child features: {0}\n".format(record))
                parents = record.attributes['Parent'].split(',') # Parent attribute could have multiple CSVs
                # expecting ordered structure within gene, although genes can be out of order across GFF file
                if self.geneFeatureArr[-1]['ID'] in parents:
                    childDict = {'ID': record.attributes['ID'],
                                 'Name': '' # display name with isoform number
                                 }
                    self.geneFeatureArr[-1]['children'].append(childDict)
                    # adding mRNA/tRNA feature record
                    self.geneFeatureArr[-1]['records'].append(record)
                else:
                    self.logger.error("gene: {0} structure unordered or unresolved parents found in gff3 record: {1}".format(self.geneFeatureArr[-1]['ID'], record))
                    raise Exception("gene: {0} structure unordered or unresolved parents found in gff3 record: {1}".format(self.geneFeatureArr[-1]['ID'], record))
            else:
                # other grand-child features
                if not 'Parent' in record.attributes:
                    self.logger.error("found gff3 feature with missing Parent attribute. Need parents for child features: {0}\n".format(record))
                    raise Exception("found gff3 feature with missing Parent attribute. Need parents for child features: {0}\n".format(record))
                parents = record.attributes['Parent'].split(',')
                if self.geneFeatureArr[-1]['children'][-1]['ID'] in parents:
                    # adding grand-child feature record
                    self.geneFeatureArr[-1]['records'].append(record)
                else:
                    self.logger.error("gene: {0} structure unordered or unresolved parents found in gff3 record: {1}".format(self.geneFeatureArr[-1]['ID'], record))
                    raise Exception("gene: {0} structure unordered or unresolved parents found in gff3 record: {1}".format(self.geneFeatureArr[-1]['ID'], record))
        
        
        # Now let's check the nomatch_seqids and assign auto numbering (starting after the max in seqidDict) & add them to seqidDict
        if self.nomatch_seqid_bin:
            # Get current max
            if self.seqidDict:
                max_seqnum = max(v for k, v in self.seqidDict.items())
            else:
                # may be nothing matched
                max_seqnum = 0
            # Add nomatch_seqids to seqidDict
            for num, nomatch_seqid in enumerate(sorted(self.nomatch_seqid_bin), start=max_seqnum+1):
                self.seqidDict[nomatch_seqid] = num
            
        
        # Update seqid in geneDict with just the number - so sorting happens correctly later
        # There are reasons this is done after parsing all records - seqidDict might be incomplete - see above step
        for geneDict in self.geneFeatureArr:
            geneDict['seqid'] = self.seqidDict[geneDict['seqid']]
                
        self.logger.info("Parsed a total of {0} genes.".format(len(self.geneFeatureArr)))
                
    
    def createDisplayName(self):
        """sort genes and create display names for gene, mRNA/tRNA features"""
        self.buildGeneSet()
        self.logger.info("Creating new display names and the ID<->Name map...")
        # sorting
        sorted_geneFeatureArr = sorted(self.geneFeatureArr, key=itemgetter('seqid', 'start'))
        # seqnum in display name
        max_seqnum = max(v for k, v in self.seqidDict.items())
        seqid_width = len(str(max_seqnum))
        # gene numbering starts from 100, 110, 120
        max_genenum = GENE_COUNTER_START + (len(sorted_geneFeatureArr) - 1)*GENE_COUNTER_STEP
        genenum_width = len(str(max_genenum))
        
        # loop through the sorted gene array
        for geneCounter, geneDict in enumerate(sorted_geneFeatureArr):
            seqnum = geneDict['seqid']
            genenum = GENE_COUNTER_START + GENE_COUNTER_STEP*geneCounter 
            # create the display name for gene feature
            # dpzm01g00100.xx
            geneDict['Name'] = "{0}{1}{2}g{3}.{4}".format(DP_PREFIX,
                                                          self.speciesCode.lower(),
                                                          str(seqnum).zfill(seqid_width),
                                                          str(genenum).zfill(genenum_width),
                                                          self.gsapRunId)
            #debug
            if geneDict['ID'] in self.id2nameMap:
                self.logger.debug("Duplicate ID: {0} feature found. Old name: {1}, New name: {2}".format(geneDict['ID'],
                                                                                             self.id2nameMap[geneDict['ID']],
                                                                                             geneDict['Name']))
            self.id2nameMap[geneDict['ID']] = geneDict['Name']
            # create the display name for child mRNA/tRNA features
            # dpzm01g00100.xx.1 - check for isoforms
            if len(geneDict['children']) == 1:
                # Only one isoform
                geneDict['children'][0]['Name'] = "{0}.{1}".format(geneDict['Name'], 1)
                #debug
                if geneDict['children'][0]['ID'] in self.id2nameMap:
                    self.logger.debug("Duplicate ID: {0} feature found. Old name: {1}, New name: {2}".format(geneDict['children'][0]['ID'],
                                                                                                 self.id2nameMap[geneDict['children'][0]['ID']],
                                                                                                 geneDict['children'][0]['Name']))
                self.id2nameMap[geneDict['children'][0]['ID']] = geneDict['children'][0]['Name']
            elif len(geneDict['children']) > 1:
                # More than one isoforms
                # sort them by ID & assign isoform number in increasing order
                sorted_children = sorted(geneDict['children'], key=lambda x: (len(x['ID']), x['ID']))
                for i, childDict in enumerate(sorted_children):
                    childDict['Name'] = "{0}.{1}".format(geneDict['Name'], i+1)
                    #debug
                    if childDict['ID'] in self.id2nameMap:
                        self.logger.debug("Duplicate ID: {0} feature found. Old name: {1}, New name: {2}".format(childDict['ID'],
                                                                                                     self.id2nameMap[childDict['ID']],
                                                                                                     childDict['Name']))
                    self.id2nameMap[childDict['ID']] = childDict['Name']
        
        self.logger.info("Created a total of {0} ID <-> display name mappings.".format(len(self.id2nameMap)))
        
    
    def writeGff3(self):
        """writes the final output sorted gff3 file with new IDs, display names and another ID<->Name mapping file"""
        self.createDisplayName()
        # writing output gff3 file
        self.logger.info("Writing gff3 records with new IDs & display names to output file: %s..." % self.outputFile)
        outputfile = open(self.outputFile, 'w')
        outputfile.write(GFF_HEADER + '\n')
        record_count = 0
        # loop through the sorted gene array - will preserve the gene structure
        for geneDict in sorted(self.geneFeatureArr, key=itemgetter('seqid', 'start')):
            for record in geneDict['records']:
                #Update source column
                record.source = SOURCE_COLUMN
                #Update display name & ID
                if record.type in ['gene', 'mRNA', 'tRNA']:
                    if 'ID' in record.attributes and record.attributes['ID'] in self.id2nameMap:
                        record.attributes['Name'] = self.id2nameMap[record.attributes['ID']]
                        record.attributes['ID'] = record.attributes['Name']
                    # debug
                    else:
                        self.logger.debug("found gff3 feature whose ID is missing or has not been mapped to name: {}\n".format(record))
                else:
                    #Update ID only - other child features
                    #format: parentId.exon1, parentId.cds1, parentId.utr5p1, parentId.utr3p1 etc
                    if 'ID' in record.attributes:
                        parts = record.attributes['ID'].split('.')
                        parentId = '.'.join(parts[:-1])
                        if parentId in self.id2nameMap:
                            record.attributes['ID'] = '.'.join([self.id2nameMap[parentId], parts[-1]])
                        # debug
                        else:
                            self.logger.debug("found gff3 feature whose parent ID {0} has not been mapped to name: {1}\n".format(parentId, record))
                #Update Parents    
                if 'Parent' in record.attributes:
                    # Parent attribute could have multiple CSVs
                    new_parents = []
                    for parent in record.attributes['Parent'].split(','):
                        if parent in self.id2nameMap:
                            new_parents.append(self.id2nameMap[parent])
                        # debug
                        else:
                            new_parents.append(parent)
                            self.logger.debug("found gff3 feature whose parent ID {0} has not been mapped to name: {1}\n".format(parent, record))
                    record.attributes['Parent'] = ','.join(new_parents)
                outputfile.write(str(record) + '\n')
                record_count += 1
            #outputfile.write('###\n')
        outputfile.close()
        self.logger.info("{0} records written to output gff3 file: {1}".format(record_count, self.outputFile))
        
        # writing seqid<->num mapping file
        #self.logger.info("Writing seqid<->num mappings to output file: %s..." % self.seqidOutputFile)
        outputfile = open(self.seqidOutputFile, 'w')
        for k, v in sorted(self.seqidDict.items(), key=itemgetter(1)):
            outputfile.write("{0}\t{1}\n".format(k, v))
        outputfile.close()
        self.logger.info("{0} seqid<-> num mappings written to output file: {1}".format(len(self.seqidDict), self.seqidOutputFile))
        
        # writing ID<->Name mapping file
        #self.logger.info("Writing ID<->Name mappings to output file: %s..." % self.id2nameOutputFile)
        outputfile = open(self.id2nameOutputFile, 'w')
        for k, v in self.id2nameMap.items():
            outputfile.write("{0}\t{1}\n".format(k, v))
        outputfile.close()
        self.logger.info("{0} ID<->Name mappings written to output file: {1}".format(len(self.id2nameMap), self.id2nameOutputFile))
        self.logger.info("Analysis all well done.")
        
    
    def speciescode_by_assemblyid(self, assemblyId):
        """ returns species code given assembly id """
        url = GAIA_ASSEMBLY_API + str(assemblyId)
        req = urllib2.Request(url)
        try:
            response = urllib2.urlopen(req)
            resultUnicode = json.loads(response.read().decode())
            result = byteify(resultUnicode)
            #print (json.dumps(resultUnicode, indent=2))
        except urllib2.HTTPError as e:
            self.logger.error("GAIA Genome Repo unreachable or couldn't fulfill the request: {0}. Error Code: {1}. Reason: {2}".format(url, e.code, e.reason))
            raise Exception("GAIA Genome Repo unreachable or couldn't fulfill the request: {0}. Error Code: {1}. Reason: {2}".format(url, e.code, e.reason))
        
        # get taxon_id and then call species code service
        if (not result) or ('taxon_id' not in result):
            self.logger.error("GAIA Genome Repo request: {0} came back empty: {1}".format(url, result))
            raise Exception("GAIA Genome Repo request: {0} came back empty: {1}".format(url, result))
        else:
            # print result['taxon_id']
            self.taxonId = result['taxon_id']
            return self.speciescode_by_taxonid(self.taxonId)

    def speciescode_by_taxonid(self, taxonId):
        """ returns species code given internal taxon id """
        url = SPECIESCODE_TAXON_API + str(taxonId)
        req = urllib2.Request(url)
        try:
            response = urllib2.urlopen(req)
            resultUnicode = json.loads(response.read().decode())
            result = byteify(resultUnicode)
            #print (json.dumps(resultUnicode, indent=2))
        except urllib2.HTTPError as e:
            self.logger.error("Species code service unreachable or couldn't fulfill the request: {0}. Error Code: {1}. Reason: {2}".format(url, e.code, e.reason))
            raise Exception("Species code service unreachable or couldn't fulfill the request: {0}. Error Code: {1}. Reason: {2}".format(url, e.code, e.reason))
        
        # get species code and return
        if (not result) or ('SpeciesCode' not in result):
            self.logger.error("Species code request: {0} came back empty: {1}".format(url, result))
            raise Exception("Species code request: {0} came back empty: {1}".format(url, result))
        elif not result['SpeciesCode']:
            self.logger.error(result['Message'])
            raise Exception(result['Message'])
        else:
            # print result['SpeciesCode']
            return result['SpeciesCode']
        
    
    def speciescode_by_sciname(self, sciname):
        """ returns species code given scientific name """
        url = urllib2.quote(SPECIESCODE_SCINAME_API + sciname, safe=":/")
        #print url
        req = urllib2.Request(url)
        try:
            response = urllib2.urlopen(req)
            resultUnicode = json.loads(response.read().decode())
            result = byteify(resultUnicode)
            #print (json.dumps(resultUnicode, indent=2))
        except urllib2.HTTPError as e:
            self.logger.error("Species code service unreachable or couldn't fulfill the request: {0}. Error Code: {1}. Reason: {2}".format(url, e.code, e.reason))
            raise Exception("Species code service unreachable or couldn't fulfill the request: {0}. Error Code: {1}. Reason: {2}".format(url, e.code, e.reason))
        
        # get species code and return
        if (not result) or ('SpeciesCode' not in result):
            self.logger.error("Species code request: {0} came back empty: {1}".format(url, result))
            raise Exception("Species code request: {0} came back empty: {1}".format(url, result))
        elif not result['SpeciesCode']:
            self.logger.error(result['Message'])
            raise Exception(result['Message'])
        else:
            # print result['SpeciesCode']
            return result['SpeciesCode']



def run():
    """run the analysis"""

    parser = argparse.ArgumentParser(
        description="Command line tool to create display names for gene and mRNA/tRNA features in GFF3",
        epilog=textwrap.dedent('''\
                                  
Example:
python %(prog)s -i input.gff3 -o /temp/output.gff3 -r 20 \
--seqid 'ZmChr(?P<num>\d+)v2' -a 219 -l /temp/log.txt \
'''),
        formatter_class=argparse.RawTextHelpFormatter)

    # Mandatory arguments
    parser.add_argument(
        '-i', dest='inputFile', required=True, help="REQUIRED: input gff3 file.")
    parser.add_argument(
        '-o', dest='outputFile', required=True, help="REQUIRED: output gff3 file. Full path.")
    parser.add_argument(
        '-r', dest='gsapRunId', required=True, help="REQUIRED: GSAP run_id.")
    parser.add_argument(
        '--seqid', dest='seqIdRegex', required=True,
        help="REQUIRED: Regex to match the seqid and extract the sequence number. \
Ex: 'ZmChr(?P<num>\d+)v2' for ZmChr01v2. \
Match the (?P<num>\d+) group, as is, with where the seq numbers appear in seqid.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-a', dest='assemblyId', type=int, help='Assembly ID from GAIA Genome Repo.')
    group.add_argument(
        '-t', dest='taxonId', type=int, help='Internal Taxon ID from Taxon Repo.')
    group.add_argument(
        '-s', dest='sciName', type=str, help='exact Scientific name. One of {-a, -t, -s} REQUIRED.')
    parser.add_argument(
        '-l', dest='logFile', help="log file, full path. Default: stderr.")

   
    args = parser.parse_args()

    # open log file first
    logFh = None
    if not args.logFile:
        logFh = sys.stderr
    else:
        logFh = open(args.logFile, 'a', 0) # no buffering.
    logger = phi.Logger.Logger(DN_ANALYSIS_NAME, logFh)

    # parsing options
    inputFile = args.inputFile
    outputFile = args.outputFile

    if outputFile[0] != '/':
        logger.error("output file %s in not in full path" % outputFile)
        sys.stderr.write("output file %s in not in full path\n" % outputFile)
        exit(1)

    
    # running command
    logger.info("The command is as below (quotes might have been removed \
when restoring the command):\n%s" % ' '.join(sys.argv))

    try:
        analysisObj = DisplayName(inputFile, outputFile, args.gsapRunId,
            args.seqIdRegex, args.assemblyId, args.taxonId, args.sciName, logger)
        analysisObj.writeGff3()

    except:
        import traceback
        stackTrace = traceback.format_exc()  # + any additional info.
        # log the error
        logger.error(stackTrace)

        raise


#! /usr/bin/env python3

# Copyright 2018, University of New Hampshire

import argparse
import gzip
import multiprocessing
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from urllib import request
from Bio import SeqIO

LOG_ERROR, LOG_WARN, LOG_INFO = range(3)
PALADIN_EXEC = "/home/anthonyw/repo_personal/unh/paladin/paladin"
RANK_LIST = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_ORDER = {rank: order for (order, rank) in enumerate(RANK_LIST)}
# All ranks used by NCBI: ["superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "superorder", "order", "parvorder", "suborder", "infraorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma"]


class FileStore:
    """ Manage temporary and cache files """
    OUTPUT_BASE = "mitobin"
    CACHE_PATH = os.path.expanduser("~/.mitobin")
    EXPIRE_AGE = 30
    FTYPE_TEMP, FTYPE_CACHE, FTYPE_OUTPUT, FTYPE_USER = range(4)
    FOPT_NORMAL, FOPT_GZIP, FOPT_TAR = range(3)

    # Entries and system generated temporary path
    entries = dict()
    temp_path = ""
    output_path = ""

    @classmethod
    def _populate(cls):
        """ Add all entries to the file store """
        FileStore("input", "input", "input.fq", None, FileStore.FTYPE_TEMP, FileStore.FOPT_NORMAL)
        FileStore("reference", "reference", "reference.faa", None, FileStore.FTYPE_CACHE, FileStore.FOPT_NORMAL)
        FileStore("sam", "sam", "{0}.sam".format(FileStore.OUTPUT_BASE), None, FileStore.FTYPE_OUTPUT, FileStore.FOPT_NORMAL)
        FileStore("synonyms", "synonyms", "synonyms.tsv", None, FileStore.FTYPE_USER, FileStore.FOPT_NORMAL)
        FileStore("lineage", "lineage", "taxdump.tar.gz", "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", FileStore.FTYPE_TEMP, FileStore.FOPT_TAR)
        FileStore("lineage-names", "lineage-names", "names.dmp", None, FileStore.FTYPE_TEMP, FileStore.FOPT_NORMAL)
        FileStore("lineage-nodes", "lineage-nodes", "nodes.dmp", None, FileStore.FTYPE_TEMP, FileStore.FOPT_NORMAL)
        FileStore("mito-gff", "mito-gff1", "mito1.gpff.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.gpff.gz", FileStore.FTYPE_TEMP, FileStore.FOPT_GZIP)
        FileStore("mito-gff", "mito-gff2", "mito2.gpff.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.protein.gpff.gz", FileStore.FTYPE_TEMP, FileStore.FOPT_GZIP)
        FileStore("mito-seq", "mito-seq1", "mito1.faa.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.faa.gz", FileStore.FTYPE_TEMP, FileStore.FOPT_GZIP)
        FileStore("mito-seq", "mito-seq2", "mito2.faa.gz", "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.protein.faa.gz", FileStore.FTYPE_TEMP, FileStore.FOPT_GZIP)

    @classmethod
    def init(cls, output_path):
        """ Initialize the file store  """
        FileStore.temp_path = tempfile.mkdtemp(prefix=FileStore.OUTPUT_BASE)
        FileStore.output_path = output_path

        # Create normal directories
        if not os.path.exists(FileStore.output_path):
            os.makedirs(FileStore.output_path)
        if not os.path.exists(FileStore.CACHE_PATH):
            os.makedirs(FileStore.CACHE_PATH)

        # Populate file store
        FileStore._populate()

    @classmethod
    def destroy(cls):
        """ Destroy the file store """
        # Close all open files
        for group in FileStore.entries.values():
            for entry in group.values():
                if entry.handle:
                    entry.handle.close()

        # Delete temporary directorty
        shutil.rmtree(FileStore.temp_path)

    @classmethod
    def get_entry(cls, name, group=None):
        """ Get a file entry from the store """
        if group:
            return FileStore.entries[group][name]
        else:
            return FileStore.entries[name][name]

    @classmethod
    def get_group(cls, group):
        """ Get all entries for a group """
        return FileStore.entries[group].values()

    def __init__(self, group, fid, name, url, ftype, options):
        self.fid = fid
        self.url = url
        self.ftype = ftype
        self.options = options
        self.handle = None

        # Set full path
        self.set_path(name, ftype)

        # Add current entry to file store
        FileStore.entries.setdefault(group, dict())
        FileStore.entries[group][fid] = self

    def set_path(self, name, ftype):
        """ Calculate and set the full path associated with the file """
        if ftype == FileStore.FTYPE_TEMP:
            self.path = os.path.join(FileStore.temp_path, name)
        if ftype == FileStore.FTYPE_CACHE:
            self.path = os.path.join(FileStore.CACHE_PATH, name)
        if ftype == FileStore.FTYPE_OUTPUT:
            self.path = os.path.join(FileStore.output_path, name)
        if ftype == FileStore.FTYPE_USER:
            self.path = name

    def get_handle(self, mode, reset=True):
        """ Return file handle and open if currently closed """
        if not self.handle:
            # Note, mode currently cannot be changed after initial open
            if self.options == FileStore.FOPT_GZIP:
                self.handle = gzip.open(self.path, mode)
            else:
                self.handle = open(self.path, mode)

        if reset:
            self.handle.seek(0)

        return self.handle

    def download(self):
        """ Download the file """
        if not self.url: return

        request.urlretrieve(self.url, self.path)

        # Extract if an archive
        if self.options == FileStore.FOPT_TAR:
            handle = tarfile.open(self.path, "r")
            handle.extractall(os.path.dirname(self.path))

    def exists(self):
        """ Check if file exists """
        return os.path.exists(self.path)

    def expired(self):
        """ Check if file is expired """
        mtime = os.path.getmtime(self.path)
        age = time.time() - mtime
        return (age / 86400) > FileStore.EXPIRE_AGE


def log(message, level):
    """ Render log messages """
    pre = ["ERROR", "WARN", "INFO"][level]
    stamp = time.strftime("%Y%m%d %H:%M:%S", time.localtime())

    print("{0} [{1}]: {2}".format(pre, stamp, message), flush=True)

    # Quit if error
    if level == 0:
        sys.exit(1)


def parse_args():
    """ Parse command line arguments """
    # Calculate default alignment options
    cores = multiprocessing.cpu_count()
    align_default = "-t {0}".format(cores)

    # Parse arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-i", dest="input", type=str, nargs="+", required=True, help="input sequence(s)")
    arg_parser.add_argument("-o", dest="output", type=str, required=True, help="output directory")
    arg_parser.add_argument("-c", dest="reference", type=str, help="custom reference")
    arg_parser.add_argument("-s", dest="sam", type=str, help="existing SAM")
    arg_parser.add_argument("-r", dest="rank", type=int, default=0, help="taxonomic rank (default: species)")
    arg_parser.add_argument("-q", dest="quality", type=int, default=0, help="mapping quality filter (default: none)")
    arg_parser.add_argument("-a", dest="align_opts", type=str, default=align_default, help="additional aligner options (default: '{0}')".format(align_default))
    arg_parser.add_argument("-f", dest="force", action="store_true", help="force rebuild of reference")

    return arg_parser.parse_args()


def prepare_workflow(args):
    """ Prepare values dependent on workflow path """
    # Existing reference skips reference generation, existing SAM skips all both
    work_ref, work_align = True, True
    if args.sam:
        work_ref, work_align = False, False
    if args.reference:
        work_ref = False

    # Directly use input file if a single input
    if len(args.input) == 1:
        FileStore.get_entry("input").set_path(args.input[0], FileStore.FTYPE_USER)

    # Redirect reference to custom path if specified
    if args.reference:
        FileStore.get_entry("reference").set_path(args.reference, FileStore.FTYPE_USER)

    # Redirect SAM to existing path if specified
    if args.sam:
        FileStore.get_entry("sam").set_path(args.sam, FileStore.FTYPE_USER)

    return work_ref, work_align


def download_dbs(force):
    """ Download sequence, annotation, and taxonomy data from NCBI """
    ref_entry = FileStore.get_entry("reference")

    if not force:
        # Check if reference exists and is still fresh
        if ref_entry.exists() and not ref_entry.expired():
            log("Fresh reference still exists, skipping preparation", LOG_INFO)
            return False

    log("Downloading databases", LOG_INFO)

    # Download every file with a URL
    for group in FileStore.entries:
        for entry in FileStore.get_group(group):
            entry.download()

    return True


def populate_synonyms():
    """ Create {synonym: official gene name} lookup, eg {CO1: COX1} """
    synonym_lookup = dict()

    for line in FileStore.get_entry("synonyms").get_handle("rt"):
        synonym, standard = line.rstrip().split("\t")
        synonym_lookup[synonym] = standard

    return synonym_lookup


def populate_accession():
    """ Create {accession: (gene name, tax ID)} lookup, eg {YP_003058231: (ND1, 123685) """
    accession_lookup = dict()

    log("Populating accession lookup", LOG_INFO)

    for entry in FileStore.get_group("mito-gff"):
        handle = entry.get_handle("rt")

        # Parse the annotations for features of interest
        for record in SeqIO.parse(handle, "genbank"):
            acc_id = record.annotations['accessions'][0]
            gene, tax_id = "unknown", None

            # Search for CDS and source features within record
            for feature in record.features:
                if feature.type == "CDS":
                    # Check for valid gene qualifier
                    if "gene" not in feature.qualifiers: continue
                    gene = feature.qualifiers["gene"][0]

                if feature.type == "source":
                    # Check for valid cross-reference qualifier
                    if "db_xref" not in feature.qualifiers: continue
                    for db in feature.qualifiers["db_xref"]:
                        if db.startswith("taxon:"):
                            tax_id = db.split(":")[1]

            # Ensure tax_id was found
            if tax_id is None: continue

            accession_lookup[acc_id] = (gene, tax_id)

    return accession_lookup


def populate_taxonomy():
    """ Create {taxid: (parent_id, rank, name)} lookup, eg {6: (335928, "genus", "Azorhizobium")} """
    taxonomy_lookup = dict()
    node_handle = FileStore.get_entry("lineage-nodes").get_handle("rt")
    name_handle = FileStore.get_entry("lineage-names").get_handle("rt")

    log("Populating taxonomy lookup", LOG_INFO)

    # Populate lineage lookup
    for node_line in node_handle:
        node_fields = [x.strip() for x in node_line.split("|")]

        # Find corresponding entry in names table
        while True:
            name_line = name_handle.readline()
            name_fields = [x.strip() for x in name_line.split("|")]

            if int(name_fields[0]) > int(node_fields[0]):
                log("Corresponding name not found for taxonomy node {0}".format(node_fields[0]), LOG_ERROR)

            if (name_fields[0] == node_fields[0]) and (name_fields[3] == "scientific name"):
                break

        taxonomy_lookup[node_fields[0]] = tuple(node_fields[1:3] + [name_fields[1].replace(" ", "_")])

    return taxonomy_lookup


def generate_lineage(tax_id, taxonomy_lookup):
    """ Generate a lineage starting at the requested taxonomy ID """
    lineage = ["unknown"] * len(RANK_LIST)
    current_id = tax_id

    # Traverse up the lineage
    while True:
        # If current ID is unknown, return lineage as is
        if current_id not in taxonomy_lookup:
            return lineage

        # Check for top of tree
        entry = taxonomy_lookup[current_id]
        if entry[0] == "1":
            break

        # Record name at supported rank
        if entry[1] in RANK_ORDER:
            lineage[RANK_ORDER[entry[1]]] = entry[2]

        # Iterate to parent
        current_id = entry[0]

    return lineage


def write_reference(accession_lookup, taxonomy_lookup, synonym_lookup):
    """ Create modified reference from original protein fasta and new headers """
    output_handle = FileStore.get_entry("reference").get_handle("wt")

    for entry in FileStore.get_group("mito-seq"):
        input_handle = entry.get_handle("rt")

        for line in input_handle:
            if line.startswith(">"):
                # Parse header
                fields = line[1:].strip().split(" ")
                acc_id, info = fields[0].split(".")[0], " ".join(fields[1:])

                # Lookup up and write new header
                gene_raw, tax_id = accession_lookup[acc_id]
                gene = synonym_lookup.get(gene_raw, "unknown")
                lineage = "|".join(generate_lineage(tax_id, taxonomy_lookup))

                # Header format "AccID|GeneID|superkingdom|kingdom|...|species supplementary"
                output_handle.write(">{0}|{1}|{2} {3}\n".format(acc_id, gene, lineage, info))
            else:
                output_handle.write(line)


def index_reference():
    """ Create aligner index for the reference """
    command = "{0} index -r 3 {1}".format(PALADIN_EXEC, FileStore.get_entry("reference").path)
    log("Indexing reference:", LOG_INFO)
    subprocess.run(command, shell=True)


def combine_inputs(input_paths):
    """ Concatenate input sequences into single file """
    # Exit if only a single file
    if len(input_paths) == 1:
        return

    # Get a handle to the temp file to be used as the concatenated input
    output_handle = FileStore.get_entry("input").get_handle("wt")

    # Iterate through all lines of each input
    for input_file in input_paths:
        file_name = os.path.split(input_file)[1]
        FileStore("inputs", file_name, input_file, None, FileStore.FTYPE_USER, FileStore.FOPT_NORMAL)
        input_entry = FileStore.get_entry(file_name, "inputs")

        if not input_entry.exists():
            log("Input file does not exist: " + input_file, LOG_ERROR)

        # Write concatenated output to temp handle
        input_handle = input_entry.get_handle("rt")
        for line in input_handle:
            output_handle.write(line)


def exec_alignment(options):
    """ Align reads to reference """
    # 1. Execute aligner for the given reads, reference, options and output
    ref_path = FileStore.get_entry("reference").path
    input_path = FileStore.get_entry("input").path
    sam_path = FileStore.get_entry("sam").path

    # Build and execute command
    log("Executing sequence alignment:", LOG_INFO)
    command = "{0} align {1} {2} {3} > {4}".format(PALADIN_EXEC, ref_path, input_path, options, sam_path)
    subprocess.run(command, shell=True)


def populate_bins(quality, rank):
    """ Filter SAM and group aligned reads into bins """
    bin_lookup = dict()
    handle = FileStore.get_entry("sam").get_handle("rt")

    for line in handle:
        sam_fields = line.split("\t")

        # Skip headers, unmapped, and low quality
        if line.startswith("@"):
            continue

        if sam_fields[1] == "4":
            continue

        if int(sam_fields[4]) < quality:
            continue

        # Confine rank to maximum if necessary
        ref_fields = sam_fields[2].split("|")
        rank_count = len(ref_fields) - 3
        if rank > rank_count:
            rank = rank_count

        # Save to lookup
        read_id = ":".join(sam_fields[0].split(":")[3:])
        bin_lookup[read_id] = (ref_fields[1], ref_fields[len(ref_fields) - rank - 1])

    return bin_lookup


def write_bins(bin_lookup):
    """ Separate and write reads into appropriate fastq file """
    read_handle = FileStore.get_entry("input").get_handle("rt")

    for line in read_handle:
        write_handle = None
        read_id = line[1:].split(" ")[0].rstrip()

        # Strip paired end strand ID
        if len(read_id) > 1 and read_id[-2] == "/":
            read_id = read_id[:-2]

        # Obtain appropriate handle and write
        if read_id in bin_lookup:
            bin_id = "{0}-{1}".format(*bin_lookup[read_id])

            # Add to FileStore if new
            if bin_id not in FileStore.entries:
                FileStore(bin_id, bin_id, "{0}.fq".format(bin_id), None, FileStore.FTYPE_OUTPUT, FileStore.FOPT_NORMAL)

            # Write header
            write_handle = FileStore.get_entry(bin_id).get_handle("wt", False)
            write_handle.write(line)

        # Write/skip subsequent 3 lines
        for _ in range(3):
            remaining = read_handle.readline()
            if write_handle:
                write_handle.write(remaining)


# Parse arguments
args = parse_args()

# Initialize the file store
FileStore.init(args.output)

# Prepare workflow steps
work_ref, work_align = prepare_workflow(args)

try:
    # Generate reference (optional)
    if work_ref:
        if download_dbs(args.force):
            accession_lookup = populate_accession()
            taxonomy_lookup = populate_taxonomy()
            synonym_lookup = populate_synonyms()

            write_reference(accession_lookup, taxonomy_lookup, synonym_lookup)
            index_reference()

    # Combine inputs (mandatory)
    combine_inputs(args.input)

    # Align reads (optional)
    if work_align:
        exec_alignment(args.align_opts)

    # Populate and write bins (mandatory)
    bin_lookup = populate_bins(args.quality, args.rank)
    write_bins(bin_lookup)

finally:
    # Destroy the file store on exit
    FileStore.destroy()

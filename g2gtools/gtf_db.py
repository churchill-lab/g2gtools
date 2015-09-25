# -*- coding: utf-8 -*-

from collections import OrderedDict
import os
import sys
try:
    from pysqlite2 import dbapi2 as sqlite3
except:
    try:
        import sqlite3
    except:
        print 'sqlite module needs to be installed'
        sys.exit()

import time

from . import G2GValueError
from .g2g_utils import configure_logging, format_time, get_logger, Location
from .gtf import GTF
import g2g_fileutils as g2g_fu

LOG = get_logger()

SQL_CREATE_GTF_TABLE = """
    CREATE TABLE gtf (
        _key INTEGER PRIMARY KEY,
        gene_id TEXT,
        transcript_id TEXT,
        ensembl_id TEXT,
        seqid TEXT NOT NULL,
        start INTEGER NOT NULL,
        end INTEGER NOT NULL,
        strand INTEGER,
        score TEXT,
        source_key INTEGER NOT NULL,
        type_key INTEGER NOT NULL,
        frame TEXT
    )
"""

SQL_CREATE_GTF_LOOKUP_TABLE = """
    CREATE TABLE gtf_lookup (
        _key INTEGER PRIMARY KEY,
        gtf_key INTEGER NOT NULL,
        attribute_key TEXT NOT NULL,
        value TEXT NOT NULL
    )
"""

SQL_CREATE_GTF_TYPES_TABLE = """
    CREATE TABLE gtf_types (
        _key INTEGER PRIMARY KEY,
        gtf_type TEXT NOT NULL
    )
"""

SQL_CREATE_GTF_SOURCES_TABLE = """
    CREATE TABLE gtf_sources (
        _key INTEGER PRIMARY KEY,
        gtf_source TEXT NOT NULL
    )
"""

SQL_CREATE_GTF_ATTRIBUTES_TABLE = """
    CREATE TABLE gtf_attributes (
        _key INTEGER PRIMARY KEY,
        gtf_attribute TEXT NOT NULL
    )
"""

SQL_INSERT_GTF_TABLE = """
    INSERT
      INTO gtf
    VALUES (null, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""

SQL_INSERT_GTF_LOOKUP_TABLE = """
    INSERT
      INTO gtf_lookup
    VALUES (null, ?, ?, ?)
"""

SQL_INSERT_GTF_TYPES_TABLE = """
    INSERT
      INTO gtf_types
    VALUES (?, ?)
"""

SQL_INSERT_GTF_SOURCES_TABLE = """
    INSERT
      INTO gtf_sources
    VALUES (?, ?)
"""

SQL_INSERT_GTF_ATTRIBUTES_TABLE = """
    INSERT
      INTO gtf_attributes
    VALUES (?, ?)
"""

SQL_INDICES_GTF = [
    "CREATE INDEX idx_gtf_gene_id ON gtf(gene_id ASC)",
    "CREATE INDEX idx_gtf_transcript_id ON gtf(transcript_id ASC)",
    "CREATE INDEX idx_gtf_ensembl_id ON gtf(ensembl_id ASC)",
    "CREATE INDEX idx_gtf_seqid ON gtf(seqid ASC)",
    "CREATE INDEX idx_gtf_start ON gtf(start ASC)",
    "CREATE INDEX idx_gtf_end ON gtf(end ASC)",
    "CREATE INDEX idx_gtf_strand ON gtf(strand ASC)",
    "CREATE INDEX idx_gtf_score ON gtf(score ASC)",
    "CREATE INDEX idx_gtf_source_key ON gtf(source_key ASC)",
    "CREATE INDEX idx_gtf_type_key ON gtf(type_key ASC)",
    "CREATE INDEX idx_gtf_frame ON gtf(frame ASC)"
]

SQL_INDICES_GTF_LOOKUP = [
    "CREATE INDEX idx_gtf_gtf_key ON gtf_lookup(gtf_key ASC)",
    "CREATE INDEX idx_gtf_lookup_attribute_key ON gtf_lookup(attribute_key ASC)",
    "CREATE INDEX idx_gtf_lookup_value ON gtf_lookup(value ASC)"
]

SQL_INDICES_GTF_TYPES = [
    "CREATE INDEX idx_gtf_types_gtf_type ON gtf_types(gtf_type ASC)",
]

SQL_INDICES_GTF_SOURCES = [
    "CREATE INDEX idx_gtf_sources_gtf_source ON gtf_sources(gtf_source ASC)",
]

SQL_INDICES_GTF_ATTRIBUTES = [
    "CREATE INDEX idx_gtf_attributes_gtf_attribute ON gtf_attributes(gtf_attribute ASC)",
]


SQL_GENE_FOR_ANY = """
SELECT g._key, g.ensembl_id, g.gene_id, g.transcript_id, g.seqid, s.gtf_source, t.gtf_type, g.seqid, g.start, g.end, g.strand, a.gtf_attribute, l.value
  FROM gtf g, gtf g2, gtf_lookup l, gtf_attributes a, gtf_sources s, gtf_types t
 WHERE g._key = l.gtf_key
   AND g.source_key = s._key
   AND g.type_key = t._key
   AND l.attribute_key = a._key
   AND g.gene_id = g2.gene_id
   and g2.ensembl_id = :ensembl_id
 ORDER BY g._key
"""

SQL_GENES_SIMPLE = """
SELECT distinct g.ensembl_id, g.seqid, g.start, g.end, g.strand
  FROM gtf g
 WHERE g.gene_id is not null
   AND g.transcript_id is null
"""
SQL_GENES_SIMPLE_ORDER_BY = " ORDER BY g._key "

# TODO need to get exons and transcripts

SQL_TRANSCRIPTS_SIMPLE = """
SELECT g._key, g.ensembl_id transcript_id, g.seqid transcript_seqid, g.start transcript_start, g.end transcript_end, g.strand transcript_strand, t.gtf_type transcript_type,
       g2._key _child_key, g2.ensembl_id child_id, g2.seqid child_seqid, g2.start child_start, g2.end child_end, g2.strand child_strand, t2.gtf_type child_type,
       g.gene_id, a.gtf_attribute, l.value
  FROM gtf g, gtf g2, gtf_types t, gtf_types t2, gtf_lookup l, gtf_attributes a
 WHERE g.type_key = t._key
   AND g.transcript_id is not null
   AND t.gtf_type = 'transcript'
   AND g2.transcript_id = g.transcript_id
   AND g2.type_key = t2._key
   AND t2.gtf_type in ('exon', 'transcript')
   AND l.gtf_key = g2._key
   AND l.attribute_key = a._key
   AND a.gtf_attribute in ('exon_number', 'gene_name')
"""

SQL_TRANSCRIPTS_SIMPLE_ORDER_BY = " ORDER BY l._key "

def count_me(conn, table):
    c = conn.cursor()
    c.execute("select count(1) from {}".format(table))
    for row in c:
        print str(row)
    c.close()

def gtf2db(input_file, output_file):
    """
    Convert a GTF file into SQLite

    :param input_file: the GTF file to convert
    :param output_file: The generated database file
    """
    start = time.time()

    input_file = g2g_fu.check_file(input_file, 'r')
    output_file = g2g_fu.check_file(output_file, 'w')

    g2g_fu.delete_file(output_file)

    LOG.info("GTF FILE: {0}".format(input_file))
    LOG.info("DB File: {0}".format(output_file))

    conn = sqlite3.connect(output_file)
    c = conn.cursor()

    LOG.debug("Generating tables")
    c.execute(SQL_CREATE_GTF_TABLE)
    c.execute(SQL_CREATE_GTF_LOOKUP_TABLE)
    c.execute(SQL_CREATE_GTF_SOURCES_TABLE)
    c.execute(SQL_CREATE_GTF_TYPES_TABLE)
    c.execute(SQL_CREATE_GTF_ATTRIBUTES_TABLE)



    gtf_types = {}
    gtf_sources = {}
    gtf_attributes = {}

    LOG.info("Parsing GTF file...")

    gtf_file = GTF(input_file)

    counter = 0

    for record in gtf_file:
        if counter and counter % 100000 == 0:
            LOG.info("Processed {0:,} records".format(counter))

        if record.type not in gtf_types:
            _type_key = len(gtf_types.keys())
            gtf_types[record.type] = _type_key
        else:
            _type_key = gtf_types[record.type]

        if record.source not in gtf_sources:
            _source_key = len(gtf_sources.keys())
            gtf_sources[record.source] = _source_key
        else:
            _source_key = gtf_sources[record.source]

        strand = 0
        if record.strand in ['+', '-']:
            strand = 1 if record.strand == '+' else -1

        gene_id = record.attributes['gene_id']
        transcript_id = record.attributes['transcript_id'] if 'transcript_id' in record.attributes else None
        ensembl_id = None

        if record.type == 'gene':
            ensembl_id = record.attributes['gene_id']
        elif record.type == 'transcript':
            ensembl_id = record.attributes['transcript_id']
        elif record.type == 'exon':
            ensembl_id = record.attributes['exon_id']
        else:
            ensembl_id = record.attributes['protein_id'] if 'protein_id' in record.attributes else None

        c.execute(SQL_INSERT_GTF_TABLE, (gene_id, transcript_id, ensembl_id, record.seqid, record.start, record.end, strand, record.score, _source_key, _type_key, record.frame))
        gtf_key = c.lastrowid

        for attribute, value in record.attributes.iteritems():
            if attribute not in ['gene_id', 'transcript_id', 'exon_id']:
                if attribute not in gtf_attributes:
                    _attribute_key = len(gtf_attributes.keys())
                    gtf_attributes[attribute] = _attribute_key
                else:
                    _attribute_key = gtf_attributes[attribute]

                c.execute(SQL_INSERT_GTF_LOOKUP_TABLE, (gtf_key, _attribute_key, value))

        counter += 1

    # save (commit) the changes
    conn.commit()

    for source, _key in gtf_sources.iteritems():
        c.execute(SQL_INSERT_GTF_SOURCES_TABLE, (_key, source))
        conn.commit()

    for type, _key in gtf_types.iteritems():
        c.execute(SQL_INSERT_GTF_TYPES_TABLE, (_key, type))
        conn.commit()

    for attribute, _key in gtf_attributes.iteritems():
        c.execute(SQL_INSERT_GTF_ATTRIBUTES_TABLE, (_key, attribute))
        conn.commit()

    LOG.info("GTF File parsed")

    LOG.info("Finalizing database...")

    for sql in SQL_INDICES_GTF:
        LOG.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_LOOKUP:
        LOG.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_TYPES:
        LOG.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_SOURCES:
        LOG.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_ATTRIBUTES:
        LOG.debug(sql)
        c.execute(sql)

    LOG.info("Database created")

    # close connection
    conn.close()

    LOG.info("Execution complete: {0}".format(format_time(start, time.time())))


class GTFObject(object):
    def __init__(self, ensembl_id=None, seqid=None, start=None, end=None, strand=None):
        self.ensembl_id = ensembl_id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.source = None
        self.name = None
        self.biotype = None


class Gene(GTFObject):
    def __init__(self, ensembl_id=None, seqid=None, start=None, end=None, strand=None):
        GTFObject.__init__(self, ensembl_id, seqid, start, end, strand)
        self.transcripts = OrderedDict()

    def __str__(self):
        return "Gene: {0} {1}:{2}-{3} ({4})".format(self.ensembl_id, self.seqid, self.start, self.end, self.strand)


class Transcript(GTFObject):
    def __init__(self, ensembl_id=None, seqid=None, start=None, end=None, strand=None):
        GTFObject.__init__(self, ensembl_id, seqid, start, end, strand)
        self.gene_ids = OrderedDict()
        self.exons = OrderedDict()

    def __str__(self):
        return "Transcript: {0} {1}:{2}-{3} ({4})".format(self.ensembl_id, self.seqid, self.start, self.end, self.strand)


class Exon(GTFObject):
    def __init__(self, ensembl_id=None, seqid=None, start=None, end=None, strand=None):
        GTFObject.__init__(self, ensembl_id, seqid, start, end, strand)
        self.gene_id = None
        self.transcript_ids = OrderedDict()
        self._exon_number = None

    @property
    def exon_number(self):
        return self._exon_number

    @exon_number.setter
    def exon_number(self, value):
        if value:
            try:
                self._exon_number = int(value)
            except ValueError, ve:
                raise G2GValueError("Illegal value for exon_number {0}, start must be an integer".format(value))

    def __str__(self):
        en = ''
        if self.exon_number:
            en = " #{0}".format(self.exon_number)
        return "Exon: {0} {1}:{2}-{3} ({4}){5}".format(self.ensembl_id, self.seqid, self.start, self.end, self.strand, en)



def location_to_sql(location, use_strand=False, overlap=True):
    sql = ""
    sql_parameters = {}

    if location:
        sql = " AND g.seqid = :seqid "
        sql_parameters = {'seqid': location.seqid}

        if overlap:
            if location.start and not location.end:
                sql = "{0} AND g.start <= :end AND g.end >= :start ".format(sql)
                sql_parameters['start'] = location.start
                sql_parameters['end'] = location.start
            elif not location.start and location.end:
                sql = "{0} AND g.start <= :end AND g.end >= :start ".format(sql)
                sql_parameters['start'] = location.end
                sql_parameters['end'] = location.end
            elif location.start and location.end:
                sql = "{0} AND g.start <= :end AND g.end >= :start ".format(sql)
                sql_parameters['start'] = location.start
                sql_parameters['end'] = location.end
        else:
            if location.start:
                sql = "{0} AND g.start >= :start ".format(sql)
                sql_parameters['start'] = location.start

            if location.end:
                sql = "{0} AND g.end <= :end ".format(sql)
                sql_parameters['end'] = location.end

        if use_strand:
            if location.strand == '-':
                strand = -1
            else:
                strand = 1

            sql = "{0} AND g.strand = :strand ".format(sql)
            sql_parameters['strand'] = strand

    return sql, sql_parameters


def get_gene(db, ensembl_id):
    """
    Retrieve a Gene object by the ensembl id.

    The ensembl id can be any value (Gene, Transcript, or Exon ID)

    """
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    cursor.execute(SQL_GENE_FOR_ANY, {'ensembl_id': ensembl_id})

    genes = OrderedDict()
    transcripts = OrderedDict()
    exons = OrderedDict()

    for r in cursor:
        _key = r['_key']
        gtf_type = r['gtf_type']
        attribute = r['gtf_attribute']
        value = r['value']

        if gtf_type == 'gene':
            gene = genes.get(_key, Gene(r['ensembl_id'], r['seqid'], r['start'], r['end'], r['strand']))

            if attribute == 'gene_name':
                gene.name = value
            elif attribute == 'gene_source':
                gene.source = value

            genes[_key] = gene

        elif gtf_type == 'transcript':
            transcript = transcripts.get(_key, Transcript(r['ensembl_id'], r['seqid'], r['start'], r['end'], r['strand']))
            transcript.gene_ids[r['gene_id']] = r['gene_id']

            if attribute == 'transcript_name':
                transcript.name = value
            elif attribute == 'transcript_source':
                transcript.source = value

            transcripts[_key] = transcript

        elif gtf_type == 'exon':
            exon = exons.get(_key, Exon(r['ensembl_id'], r['seqid'], r['start'], r['end'], r['strand']))
            exon.gene_ids = r['gene_id']
            exon.transcript_ids[r['transcript_id']] = r['transcript_id']

            if attribute == 'exon_number':
                exon.exon_number = value

            exons[_key] = exon

    genes = {gene.ensembl_id: gene for i, gene in genes.iteritems()}
    transcripts = {transcript.ensembl_id: transcript for i, transcript in transcripts.iteritems()}

    for _id, exon in exons.iteritems():
        for _tid in exon.transcript_ids:
            transcripts[_tid].exons[exon.ensembl_id] = exon

    for _id, transcript in transcripts.iteritems():
        for _gid in transcript.gene_ids:
            genes[_gid].transcripts[_id] = transcript

    cursor.close()
    conn.close()
    return genes


def get_genes_ids(db, location=None, use_strand=False, overlap=True):
    """
    Get ensembl ids for genes
    """

    sql, sql_parameters = location_to_sql(location, use_strand, overlap)
    sql = "{0} {1} {2}".format(SQL_GENES_SIMPLE, sql, SQL_GENES_SIMPLE_ORDER_BY)

    LOG.debug("SQL:\n{0}".format(sql))
    LOG.debug("PARAMETERS: {0}".format(sql_parameters))

    conn = sqlite3.connect(db)
    sqlite3.enable_callback_tracebacks(True)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    cursor.execute(sql, sql_parameters)

    gene_ids = []
    for r in cursor:
        gene_ids.append(r['ensembl_id'])

    cursor.close()
    conn.close()
    return gene_ids



def get_genes(db, location=None, use_strand=False, overlap=True):
    """
    Get Gene objects
    """
    gene_ids = get_genes_ids(db, location, use_strand, overlap)

    genes = {}
    for i, gene_id in enumerate(gene_ids):
        genes_temp = get_gene(db, gene_id)
        if genes_temp:
            for ensembl_id, gene in genes_temp.iteritems():
                genes[ensembl_id] = gene

    return genes


def get_genes_simple(db, location=None, use_strand=False, overlap=True):

    sql, sql_parameters = location_to_sql(location, use_strand, overlap)
    sql = "{0} {1} {2}".format(SQL_GENES_SIMPLE, sql, SQL_GENES_SIMPLE_ORDER_BY)

    LOG.debug("SQL:\n{0}".format(sql))
    LOG.debug("PARAMETERS: {0}".format(sql_parameters))

    conn = sqlite3.connect(db)
    sqlite3.enable_callback_tracebacks(True)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    cursor.execute(sql, sql_parameters)

    genes = []
    for r in cursor:
        genes.append(Gene(r['ensembl_id'], r['seqid'], r['start'], r['end'], r['strand']))

    cursor.close()
    conn.close()

    return genes

def get_transcripts_simple(db, location=None, use_strand=False, overlap=True):
    sql, sql_parameters = location_to_sql(location, use_strand, overlap)
    sql = "{0} {1} {2}".format(SQL_TRANSCRIPTS_SIMPLE, sql, SQL_TRANSCRIPTS_SIMPLE_ORDER_BY)

    LOG.debug("SQL:\n{0}".format(sql))
    LOG.debug("PARAMETERS: {0}".format(sql_parameters))

    conn = sqlite3.connect(db)
    sqlite3.enable_callback_tracebacks(True)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    cursor.execute(sql, sql_parameters)

    transcripts = OrderedDict()
    exons = OrderedDict()
    for r in cursor:

        if r['transcript_id'] == r['child_id']:
            # transcript
            if r['_key'] not in transcripts:
                transcripts[r['_key']] = Transcript(r['transcript_id'], r['transcript_seqid'], r['transcript_start'], r['transcript_end'], r['transcript_strand'])
        else:
            # exon
            exon = exons.get(r['_child_key'], Exon(r['child_id'], r['child_seqid'], r['child_start'], r['child_end'], r['child_strand']))
            exon.gene_id = r['transcript_id']
            exon.transcript_ids[r['transcript_id']] = r['transcript_id']
            attribute = r['gtf_attribute']
            value = r['value']

            if attribute == 'exon_number':
                exon.exon_number = value

            exons[r['_child_key']] = exon

    transcripts = {transcript.ensembl_id: transcript for i, transcript in transcripts.iteritems()}

    for _id, exon in exons.iteritems():
        for _tid in exon.transcript_ids:
            transcripts[_tid].exons[exon.ensembl_id] = exon

    cursor.close()
    conn.close()

    return transcripts.values()


if __name__ == '__main__':
    import sys
    configure_logging(5)

    genes = get_gene(sys.argv[1], sys.argv[2])
    for k,v in genes.iteritems():
        print v
        for k1, v1 in v.transcripts.iteritems():
            print "\t", v1
            for k2,v2 in v1.exons.iteritems():
                print "\t\t", v2

    print '*'*80

    genes = get_genes(sys.argv[1], Location('1', 21250000, 21500000, '-'), overlap=True, use_strand=True)
    for k,v in genes.iteritems():
        print v
        for k1, v1 in v.transcripts.iteritems():
            print "\t", v1
            for k2,v2 in v1.exons.iteritems():
                print "\t\t", v2

    print len(genes)

    """
    gene_ids = get_genes_ids(sys.argv[1], Location('1', 21250000, 21500000), overlap=True)
    print len(gene_ids)

    genes = get_genes_simple(sys.argv[1], Location('1', 21250000, 21500000, '-'), overlap=False, use_strand=True)

    for i, gene in enumerate(genes):
        print gene

    for i, transcript in enumerate(get_transcripts_simple(sys.argv[1], Location('12', 4833175, 4866873), overlap=False)):
            print transcript, len(transcript.exons)
            for k2,v2 in transcript.exons.iteritems():
                print "\t", v2

    """
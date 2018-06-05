# -*- coding: utf-8 -*-

from collections import OrderedDict
from future.utils import viewitems
import sys
import time


try:
    import sqlite3
except:
    print('sqlite module needs to be installed')
    sys.exit()

from . import exceptions
from . import g2g
from . import g2g_utils
from . import gtf

LOG = g2g.get_logger()

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
SELECT  *,
        (SELECT value 
           FROM gtf_lookup l,
                gtf_attributes a 
          WHERE l.gtf_key = g._key 
            AND l.attribute_key = a._key
            AND a.gtf_attribute = 'exon_number') AS exon_number
  FROM gtf g,
       gtf_types t
 WHERE g.type_key = t._key
   AND t.gtf_type in ('exon', 'transcript', 'gene')
"""

SQL_TRANSCRIPTS_SIMPLE_ORDER_BY = " ORDER BY g._key "


def dictify_row(cursor, row):
    """Turns the given row into a dictionary where the keys are the column names.

    Args:
        cursor (sqlite3.Cursor): the database cursor
        row (dict): the current row

    Returns:
        OrderedDict: dictionary with keys being column names
    """
    d = OrderedDict()
    for i, col in enumerate(cursor.description):
        d[col[0]] = row[i]
    return d


def dictify_cursor(cursor):
    """
    Converts all cursor rows into dictionaries where keys are the column names.

    Args:
        cursor (sqlite3.Cursor): the database cursor

    Returns:
        list: list of ``OrderedDict`` objects with keys being column names
    """
    return [dictify_row(cursor, row) for row in cursor]


def count_me(conn, table):
    c = conn.cursor()
    c.execute("select count(1) from {}".format(table))
    for row in c:
        print(str(row))
    c.close()


def gtf2db(input_file, output_file):
    """
    Convert a GTF file into SQLite

    :param input_file: the GTF file to convert
    :param output_file: The generated database file
    """
    start = time.time()

    input_file = g2g_utils.check_file(input_file, 'r')
    output_file = g2g_utils.check_file(output_file, 'w')

    g2g_utils.delete_file(output_file)

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

    gtf_file = gtf.GTF(input_file)

    counter = 0
    prev_gene_id = None

    for record in gtf_file:
        LOG.debug('LINE={}'.format(record))
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

        if prev_gene_id != gene_id:
            LOG.debug('-'*80)
            prev_gene_id = gene_id

        transcript_id = record.attributes['transcript_id'] if 'transcript_id' in record.attributes else None
        ensembl_id = None

        LOG.debug('transcript_id = {}'.format(transcript_id))

        if record.type == 'gene':
            ensembl_id = record.attributes['gene_id']
        elif record.type == 'transcript':
            ensembl_id = record.attributes['transcript_id']
        elif record.type == 'exon':
            ensembl_id = record.attributes['exon_id']
        else:
            ensembl_id = record.attributes['protein_id'] if 'protein_id' in record.attributes else None

        LOG.debug('ensembl_id = {}'.format(ensembl_id))

        c.execute(SQL_INSERT_GTF_TABLE, (gene_id, transcript_id, ensembl_id, record.seqid, record.start, record.end, strand, record.score, _source_key, _type_key, record.frame))

        LOG.debug('INSERTING GTF = {}'.format(str((gene_id, transcript_id, ensembl_id, record.seqid, record.start, record.end, strand, record.score, _source_key, _type_key, record.frame))))
        gtf_key = c.lastrowid
        LOG.debug('gtf_key={}'.format(gtf_key))

        for (attribute, value) in viewitems(record.attributes):
            if attribute not in ['gene_id', 'transcript_id', 'exon_id']:
                if attribute not in gtf_attributes:
                    _attribute_key = len(gtf_attributes.keys())
                    gtf_attributes[attribute] = _attribute_key
                else:
                    _attribute_key = gtf_attributes[attribute]

                LOG.debug('inserting {}'.format(str((gtf_key, _attribute_key, value))))
                c.execute(SQL_INSERT_GTF_LOOKUP_TABLE, (gtf_key, _attribute_key, value))

        counter += 1

    # save (commit) the changes
    conn.commit()

    for (source, _key) in viewitems(gtf_sources):
        c.execute(SQL_INSERT_GTF_SOURCES_TABLE, (_key, source))
        conn.commit()

    for (type, _key) in viewitems(gtf_types):
        c.execute(SQL_INSERT_GTF_TYPES_TABLE, (_key, type))
        conn.commit()

    for (attribute, _key) in viewitems(gtf_attributes):
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

    LOG.info("Execution complete: {0}".format(g2g_utils.format_time(start, time.time())))


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
            except ValueError as ve:
                raise exceptions.G2GValueError("Illegal value for exon_number {0}, start must be an integer".format(value))

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
    conn.text_factory = str
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

    genes = {gene.ensembl_id: gene for i, gene in genes.items()}
    transcripts = {transcript.ensembl_id: transcript for i, transcript in transcripts.items()}

    for _id, exon in exons.items():
        for _tid in exon.transcript_ids:
            transcripts[_tid].exons[exon.ensembl_id] = exon

    for _id, transcript in transcripts.items():
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
    conn.text_factory = str
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
            for ensembl_id, gene in genes_temp.items():
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
    conn.text_factory = str
    cursor = conn.cursor()

    cursor.execute(sql, sql_parameters)

    genes = []
    for r in cursor:
        genes.append(Gene(r['ensembl_id'], r['seqid'], r['start'], r['end'], r['strand']))

    cursor.close()
    conn.close()

    return genes

#def get_transcripts_simple(db, location=None, use_strand=False, overlap=True):
def get_transcripts_simple(db):
    sql = "{0} {1}".format(SQL_TRANSCRIPTS_SIMPLE, SQL_TRANSCRIPTS_SIMPLE_ORDER_BY)

    LOG.debug("SQL:\n{0}".format(sql))

    conn = sqlite3.connect(db)

    #LOG.debug('SQLite Version: {}'.format(sqlite3.sqlite_version))
    #LOG.debug('sqlite3 driver {}, version: {}'.format(sqlite3.__name__, sqlite3.version))
    sqlite3.enable_callback_tracebacks(True)
    conn.row_factory = sqlite3.Row
    conn.text_factory = str
    cursor = conn.cursor()
    cursor.execute(sql)

    transcripts = OrderedDict()

    counter = 0
    for r in cursor:
        counter += 1
        if counter % 10000 == 0:
            LOG.debug("Parsed {} elements".format(counter))

        if r['transcript_id'] == r['ensembl_id']:
            # transcript
            if r['transcript_id'] not in transcripts:
                #LOG.debug('adding transcript {}'.format(r['transcript_id']))

                # ** deactivated
                #transcripts[r['ensembl_id']] = Transcript(r['ensembl_id'], r['seqid'], r['start'], r['end'], r['strand'])

                ## ** New method added : Since the condition i.e r['transcript_id'] == r['ensembl_id'] is true ..
                # .. and also that we are building transcript-ids, using 'transcript_id' as key is more comprehensive.
                # Also, it would be comprehensive to have : Transcript(r['transcript_id'], ......
                transcripts[r['transcript_id']] = Transcript(r['ensembl_id'], r['seqid'], r['start'], r['end'],
                                                          r['strand'])

                #LOG.debug(transcripts[r['ensembl_id']])

        # ** deactivated
        #elif r['transcript_id'] is not None and (r['gene_id'] != r['ensembl_id']):

        # ** added to prevent inserttion of 'None' key in exon.ensembl_id ??
        elif r['transcript_id'] is not None and r['ensembl_id'] is not None \
                    and (r['gene_id'] != r['ensembl_id']):
            # exon
            LOG.debug("ADDING exon")
            LOG.debug('{}:{}'.format(r['ensembl_id'], r['exon_number']))
            exon = Exon(r['ensembl_id'], r['seqid'], r['start'], r['end'], r['strand'])
            exon.gene_id = r['gene_id']
            exon.transcript_ids[r['transcript_id']] = r['transcript_id']
            exon.exon_number = r['exon_number']


            transcripts[r['transcript_id']].exons[r['ensembl_id']] =  exon


            ## ** Addition
            # When there are sinlge Exon genes we are missing transcript_id based on previous two condition

            ## A small snippet from gtftodb (sqlite file) for the gene that has only one exon
            # gene_id       transcript_id       ensembl_id
            # "AL1G10030_L"   null    "AL1G10030_L"
            # "AL1G10030_R"   null    "AL1G10030_R"
            # "AL1G10030_L"    "AL1G10030.t1_L"    "AL1G10030.t1.exon1_L"
            # "AL1G10030_R"    "AL1G10030.t1_R"    "AL1G10030.t1.exon1_R"
            # "AL1G10030_L"    "AL1G10030.t1_L"   null
            # "AL1G10030_R"    "AL1G10030.t1_R"   null
            # "AL1G10030_L"    "AL1G10030.t1_L"   null
            # "AL1G10030_R"    "AL1G10030.t1_R"   null
            # "AL1G10030_L"    "AL1G10030.t1_L"   null
            # "AL1G10030_R"    "AL1G10030.t1_R"   null

            # to include missing transcript_ids database from single exon genes adding another condition
            # transcript
            if r['transcript_id'] not in transcripts:
                # LOG.debug('adding transcript {}'.format(r['transcript_id']))

                ## ** New method: Since the condition i.e r['transcript_id'] != r['ensembl_id'] ..
                # .. and also that we are building transcript-ids, using 'transcript_id' as key is more comprehensive
                transcripts[r['transcript_id']] = Transcript(r['transcript_id'], r['seqid'], r['start'], r['end'],
                                                             r['strand'])

        else:
            #LOG.debug('gene')
            pass


    LOG.debug("Simplifying transcripts")
    #transcripts = {transcript.ensembl_id: transcript for i, transcript in transcripts.items()}
    #LOG.debug(transcripts)
    for i, transcript in transcripts.items():
        LOG.debug("Transcript={0}".format(transcript))

        for ensembl_id, exon in transcript.exons.items():
            LOG.debug("Exon ID={0};{1}".format(ensembl_id, exon))

    #for _id, exon in exons.items():
    #    LOG.debug('_id={}\texon={}'.format(_id, str(exon)))
    #    for _tid in exon.transcript_ids:
    #        #LOG.debug('_tid={}'.format(_id))
    #        LOG.debug(transcripts[_tid])
    #        transcripts[_tid].exons[exon.ensembl_id] = exon

    cursor.close()
    conn.close()

    LOG.debug("Number of transcripts: {}".format(len(transcripts.values())))

    return transcripts.values()

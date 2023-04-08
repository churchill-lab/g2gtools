# standard library imports
from collections import OrderedDict
import sqlite3
import time

# 3rd party library imports
# none

# local library imports
from .exceptions import G2GValueError
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
    "CREATE INDEX idx_gtf_frame ON gtf(frame ASC)",
]

SQL_INDICES_GTF_LOOKUP = [
    "CREATE INDEX idx_gtf_gtf_key ON gtf_lookup(gtf_key ASC)",
    "CREATE INDEX idx_gtf_lookup_attr_key ON gtf_lookup(attribute_key ASC)",
    "CREATE INDEX idx_gtf_lookup_value ON gtf_lookup(value ASC)",
]

SQL_INDICES_GTF_TYPES = [
    "CREATE INDEX idx_gtf_types_gtf_type ON gtf_types(gtf_type ASC)",
]

SQL_INDICES_GTF_SOURCES = [
    "CREATE INDEX idx_gtf_sources_gtf_source ON gtf_sources(gtf_source ASC)",
]

SQL_INDICES_GTF_ATTRIBUTES = [
    "CREATE INDEX idx_gtf_attr_gtf_att ON gtf_attributes(gtf_attribute ASC)",
]

SQL_GENE_FOR_ANY = """
SELECT g._key, g.ensembl_id, g.gene_id, g.transcript_id, g.seqid, 
       s.gtf_source, t.gtf_type, g.seqid, g.start, g.end, g.strand, 
       a.gtf_attribute, l.value
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


def gtf2db(
    gtf_file_name: str, database_file_name: str, debug_level: int | None = 0
) -> None:
    """
    Convert a GTF file into SQLite

    Args:
        gtf_file_name: The GTF file to parse.
        database_file_name: The SQLite file name to create.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    start = time.time()
    logger = g2g.get_logger(debug_level)

    gtf_file_name = g2g_utils.check_file(gtf_file_name, "r")
    database_file_name = g2g_utils.check_file(database_file_name, "w")

    g2g_utils.delete_file(database_file_name)

    logger.warn(f"GTF FILE: {gtf_file_name}")
    logger.warn(f"DB File: {database_file_name}")

    conn = sqlite3.connect(database_file_name)
    c = conn.cursor()

    logger.warn("Generating tables...")
    c.execute(SQL_CREATE_GTF_TABLE)
    c.execute(SQL_CREATE_GTF_LOOKUP_TABLE)
    c.execute(SQL_CREATE_GTF_SOURCES_TABLE)
    c.execute(SQL_CREATE_GTF_TYPES_TABLE)
    c.execute(SQL_CREATE_GTF_ATTRIBUTES_TABLE)

    gtf_types = {}
    gtf_sources = {}
    gtf_attributes = {}

    logger.warn("Parsing GTF file...")

    gtf_file = gtf.GTF(gtf_file_name)

    counter = 0
    prev_gene_id = None

    for record in gtf_file:
        logger.debug(f"LINE={record}")

        if counter and counter % 100000 == 0:
            logger.info(f"Processed {counter:,} records")

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
        if record.strand in ["+", "-"]:
            strand = 1 if record.strand == "+" else -1

        gene_id = record.attributes["gene_id"]

        if prev_gene_id != gene_id:
            logger.debug("-" * 80)
            prev_gene_id = gene_id

        if "transcript_id" in record.attributes:
            transcript_id = record.attributes["transcript_id"]
        else:
            transcript_id = None

        logger.debug(f"transcript_id = {transcript_id}")

        if record.type == "gene":
            ensembl_id = record.attributes["gene_id"]
        elif record.type == "transcript":
            ensembl_id = record.attributes["transcript_id"]
        elif record.type == "exon":
            ensembl_id = record.attributes["exon_id"]
        else:
            if "protein_id" in record.attributes:
                ensembl_id = record.attributes["protein_id"]
            else:
                ensembl_id = None

        logger.debug(f"ensembl_id = {ensembl_id}")

        c.execute(
            SQL_INSERT_GTF_TABLE,
            (
                gene_id,
                transcript_id,
                ensembl_id,
                record.seqid,
                record.start,
                record.end,
                strand,
                record.score,
                _source_key,
                _type_key,
                record.frame,
            ),
        )

        gtf_key = c.lastrowid
        logger.debug(f"gtf_key={gtf_key}")

        for attribute, value in record.attributes.items():
            if attribute not in ["gene_id", "transcript_id", "exon_id"]:
                if attribute not in gtf_attributes:
                    _attribute_key = len(gtf_attributes.keys())
                    gtf_attributes[attribute] = _attribute_key
                else:
                    _attribute_key = gtf_attributes[attribute]

                logger.debug(f"inserting {gtf_key}, {_attribute_key}, {value}")
                c.execute(
                    SQL_INSERT_GTF_LOOKUP_TABLE,
                    (gtf_key, _attribute_key, value)
                )

        counter += 1

    # save (commit) the changes
    conn.commit()

    for source, _key in gtf_sources.items():
        c.execute(SQL_INSERT_GTF_SOURCES_TABLE, (_key, source))
        conn.commit()

    for _type, _key in gtf_types.items():
        c.execute(SQL_INSERT_GTF_TYPES_TABLE, (_key, _type))
        conn.commit()

    for attribute, _key in gtf_attributes.items():
        c.execute(SQL_INSERT_GTF_ATTRIBUTES_TABLE, (_key, attribute))
        conn.commit()

    logger.warn("GTF File parsed")

    logger.warn("Finalizing database...")

    for sql in SQL_INDICES_GTF:
        logger.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_LOOKUP:
        logger.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_TYPES:
        logger.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_SOURCES:
        logger.debug(sql)
        c.execute(sql)

    for sql in SQL_INDICES_GTF_ATTRIBUTES:
        logger.debug(sql)
        c.execute(sql)

    logger.warn("Database created")

    # close connection
    conn.close()

    fmt_time = g2g_utils.format_time(start, time.time())
    logger.warn(f"Execution complete: {fmt_time}")


class GTFObject(object):
    def __init__(
        self,
        ensembl_id: str | None = None,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
        strand: str | None = None,
    ):
        self.ensembl_id = ensembl_id
        self.seqid = seqid
        self.start = start
        self.end = end
        self.strand = strand
        self.source = None
        self.name = None
        self.biotype = None


class Gene(GTFObject):
    def __init__(
        self,
        ensembl_id: str | None = None,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
        strand: str | None = None,
    ):
        GTFObject.__init__(self, ensembl_id, seqid, start, end, strand)
        self.transcripts = OrderedDict()

    def __str__(self):
        return (
            f"Gene: {self.ensembl_id} "
            f"{self.seqid}:{self.start}-{self.end} ({self.strand})"
        )


class Transcript(GTFObject):
    def __init__(
        self,
        ensembl_id: str | None = None,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
        strand: str | None = None,
    ):
        GTFObject.__init__(self, ensembl_id, seqid, start, end, strand)
        self.gene_ids = OrderedDict()
        self.exons = OrderedDict()

    def __str__(self):
        return (
            f"Transcript: {self.ensembl_id} "
            f"{self.seqid}:{self.start}-{self.end} ({self.strand})"
        )


class Exon(GTFObject):
    def __init__(
        self,
        ensembl_id: str | None = None,
        seqid: str | None = None,
        start: int | None = None,
        end: int | None = None,
        strand: str | None = None,
    ):
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
            except ValueError:
                raise G2GValueError(
                    f"Illegal value for exon_number {value}, must be an integer"
                )

    def __str__(self):
        en = ""

        if self.exon_number:
            en = " #{0}".format(self.exon_number)

        return (
            f"Exon: {self.ensembl_id} "
            f"{self.seqid}:{self.start}-{self.end} ({self.strand}) {en}"
        )


def location_to_sql(
    location: g2g.Region, use_strand: bool = False, overlap: bool | None = True
) -> tuple[str, dict]:
    """
    Utility function to convert a Region into a SQL condition.

    Args:
        location: The region.
        use_strand: Include where condition for strand.
        overlap: Include where condition for overlap.

    Returns:
        A SQL condition for the WHERE clause and the associated parameters.
    """
    sql = ""
    sql_parameters = {}

    if location:
        sql = " AND g.seqid = :seqid "
        sql_parameters = {"seqid": location.seq_id}

        if overlap:
            if location.start and not location.end:
                sql = f"{sql} AND g.start <= :end AND g.end >= :start "
                sql_parameters["start"] = location.start
                sql_parameters["end"] = location.start
            elif not location.start and location.end:
                sql = f"{sql} AND g.start <= :end AND g.end >= :start "
                sql_parameters["start"] = location.end
                sql_parameters["end"] = location.end
            elif location.start and location.end:
                sql = f"{sql} AND g.start <= :end AND g.end >= :start "
                sql_parameters["start"] = location.start
                sql_parameters["end"] = location.end
        else:
            if location.start:
                sql = f"{sql} AND g.start >= :start "
                sql_parameters["start"] = location.start

            if location.end:
                sql = f"{sql} AND g.end <= :end "
                sql_parameters["end"] = location.end

        if use_strand:
            if location.strand == "-":
                strand = -1
            else:
                strand = 1

            sql = f"{sql} AND g.strand = :strand "
            sql_parameters["strand"] = strand

    return sql, sql_parameters


def get_gene(database_file_name: str, ensembl_id: str) -> dict[str, Gene]:
    """
    Retrieve a Gene object by the ensembl id.

    The ensembl id can be any value (Gene, Transcript, or Exon ID)

    Args:
        database_file_name: Name of the SQLite database to use.
        ensembl_id: The Ensembl identifier.

    Returns:
        A dictionary of matching ids of Gene objects.
    """
    conn = sqlite3.connect(database_file_name)
    conn.row_factory = sqlite3.Row
    conn.text_factory = str
    cursor = conn.cursor()
    cursor.execute(SQL_GENE_FOR_ANY, {"ensembl_id": ensembl_id})

    genes = OrderedDict()
    transcripts = OrderedDict()
    exons = OrderedDict()

    for r in cursor:
        _key = r["_key"]
        gtf_type = r["gtf_type"]
        attribute = r["gtf_attribute"]
        value = r["value"]

        if gtf_type == "gene":
            gene = genes.get(
                _key,
                Gene(
                    r["ensembl_id"],
                    r["seqid"],
                    r["start"],
                    r["end"],
                    r["strand"]
                ),
            )

            if attribute == "gene_name":
                gene.name = value
            elif attribute == "gene_source":
                gene.source = value

            genes[_key] = gene

        elif gtf_type == "transcript":
            transcript = transcripts.get(
                _key,
                Transcript(
                    r["ensembl_id"],
                    r["seqid"],
                    r["start"],
                    r["end"],
                    r["strand"]
                ),
            )
            transcript.gene_ids[r["gene_id"]] = r["gene_id"]

            if attribute == "transcript_name":
                transcript.name = value
            elif attribute == "transcript_source":
                transcript.source = value

            transcripts[_key] = transcript

        elif gtf_type == "exon":
            exon = exons.get(
                _key,
                Exon(
                    r["ensembl_id"],
                    r["seqid"],
                    r["start"],
                    r["end"],
                    r["strand"]
                ),
            )
            exon.gene_ids = r["gene_id"]
            exon.transcript_ids[r["transcript_id"]] = r["transcript_id"]

            if attribute == "exon_number":
                exon.exon_number = value

            exons[_key] = exon

    genes = {gene.ensembl_id: gene for i, gene in genes.items()}
    transcripts = {
        tr.ensembl_id: tr for i, tr in transcripts.items()
    }

    for _id, exon in exons.items():
        for _tid in exon.transcript_ids:
            transcripts[_tid].exons[exon.ensembl_id] = exon

    for _id, transcript in transcripts.items():
        for _gid in transcript.gene_ids:
            genes[_gid].transcripts[_id] = transcript

    cursor.close()
    conn.close()

    return genes


def get_genes_ids(
    database_file_name: str,
    location: g2g.Region,
    use_strand: bool | None = False,
    overlap: bool | None = True,
) -> list[str]:
    """
    Get ensembl ids for genes.

    Args:
        database_file_name: Name of the SQLite database to use.
        location: Location to restrict to.
        use_strand: True to use strand, false to not.
        overlap: True to get genes that overlap, False otherwise.

    Returns:
        A list of Ensembl gene identifiers.
    """
    sql, sql_parameters = location_to_sql(location, use_strand, overlap)
    sql = f"{SQL_GENES_SIMPLE} {sql} {SQL_GENES_SIMPLE_ORDER_BY}"

    conn = sqlite3.connect(database_file_name)
    sqlite3.enable_callback_tracebacks(True)
    conn.row_factory = sqlite3.Row
    conn.text_factory = str
    cursor = conn.cursor()

    cursor.execute(sql, sql_parameters)

    gene_ids: list[str] = []

    for r in cursor:
        gene_ids.append(r["ensembl_id"])

    cursor.close()
    conn.close()

    return gene_ids


def get_genes(
    database_file_name: str,
    location: g2g.Region,
    use_strand: bool | None = False,
    overlap: bool | None = True,
) -> dict[str, Gene]:
    """
    Get Gene objects.

    Args:
        database_file_name: Name of the SQLite database to use.
        location: Location to restrict to.
        use_strand: True to use strand, false to not.
        overlap: True to get genes that overlap, False otherwise.
    """
    gene_ids = get_genes_ids(database_file_name, location, use_strand, overlap)

    genes: dict[str, Gene] = {}
    for i, gene_id in enumerate(gene_ids):
        genes_temp = get_gene(database_file_name, gene_id)
        if genes_temp:
            for ensembl_id, gene in genes_temp.items():
                genes[ensembl_id] = gene

    return genes


def get_genes_simple(
    database_file_name: str,
    location: g2g.Region | None = None,
    use_strand: bool | None = False,
    overlap: bool | None = True,
    debug_level: int | None = 0,
) -> list[Gene]:
    """
    Get Gene objects.

    Args:
        database_file_name: Name of the SQLite database to use.
        location: Location to restrict to.
        use_strand: True to use strand, false to not.
        overlap: True to get genes that overlap, False otherwise.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).

    Returns:
          A list of matching Gene objects.
    """
    sql, sql_parameters = location_to_sql(location, use_strand, overlap)
    sql = f"{SQL_GENES_SIMPLE} {sql} {SQL_GENES_SIMPLE_ORDER_BY}"

    logger = g2g.get_logger(debug_level)
    logger.debug(f"SQL:\n{sql}")
    logger.debug(f"PARAMETERS: {sql_parameters}")

    conn = sqlite3.connect(database_file_name)
    sqlite3.enable_callback_tracebacks(True)
    conn.row_factory = sqlite3.Row
    conn.text_factory = str
    cursor = conn.cursor()

    cursor.execute(sql, sql_parameters)

    genes: list[Gene] = []
    for r in cursor:
        genes.append(
            Gene(r["ensembl_id"], r["seqid"], r["start"], r["end"], r["strand"])
        )

    cursor.close()
    conn.close()

    return genes


def get_transcripts_simple(
    database_file_name: str, debug_level: int | None = 0
) -> list[Transcript]:
    """
    Get all the transcripts from the database.

    Args:
        database_file_name: Name of the SQLite database to use.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).

    Returns:
        A list of Transcript objects.
    """
    sql = f"{SQL_TRANSCRIPTS_SIMPLE} {SQL_TRANSCRIPTS_SIMPLE_ORDER_BY}"

    logger = g2g.get_logger(debug_level)
    logger.debug(f"SQL:\n{sql}")

    # logger.debug(f"SQLite Version: {sqlite3.sqlite_version}")
    # logger.debug(f"SQLite Driver: {sqlite3.__name__}")
    # logger.debug(f"SQLite Driver Version: {sqlite3.version}")

    conn = sqlite3.connect(database_file_name)
    sqlite3.enable_callback_tracebacks(True)
    conn.row_factory = sqlite3.Row
    conn.text_factory = str
    cursor = conn.cursor()
    cursor.execute(sql)

    transcripts: dict[str, Transcript] = OrderedDict()

    counter = 0
    for r in cursor:
        counter += 1
        if counter % 10000 == 0:
            logger.debug(f"Parsed {counter} elements")

        if r["transcript_id"] == r["ensembl_id"]:
            # transcript
            if r["transcript_id"] not in transcripts:
                # logger.debug(f"Adding transcript {r["transcript_id"]}")

                transcripts[r["transcript_id"]] = Transcript(
                    r["ensembl_id"],
                    r["seqid"],
                    r["start"],
                    r["end"],
                    r["strand"]
                )

        elif (
            r["transcript_id"] is not None
            and r["ensembl_id"] is not None
            and (r["gene_id"] != r["ensembl_id"])
        ):
            # exon
            # logger.debug("ADDING exon")
            # logger.debug(f"{r['ensembl_id']}:{r['exon_number']}")

            exon = Exon(
                r["ensembl_id"],
                r["seqid"],
                r["start"],
                r["end"],
                r["strand"]
            )

            exon.gene_id = r["gene_id"]
            exon.transcript_ids[r["transcript_id"]] = r["transcript_id"]
            exon.exon_number = r["exon_number"]

            transcripts[r["transcript_id"]].exons[r["ensembl_id"]] = exon

            # A small snippet from gtf2db (sqlite file) for the gene
            # that has only one exon
            #
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

            # to include missing transcript_ids database from single exon genes
            # adding another condition transcript
            if r["transcript_id"] not in transcripts:
                # logger.debug(f"Adding transcript {r["transcript_id"]}"))

                transcripts[r["transcript_id"]] = Transcript(
                    r["transcript_id"],
                    r["seqid"],
                    r["start"],
                    r["end"],
                    r["strand"]
                )
        else:
            # logger.debug("gene")
            pass

    logger.debug("Simplifying transcripts")

    for i, transcript in transcripts.items():
        logger.debug(f"Transcript={transcript}")

        for ensembl_id, exon in transcript.exons.items():
            logger.debug(f"Exon ID={ensembl_id};{exon}")

    cursor.close()
    conn.close()

    logger.debug(f"Number of transcripts: {len(transcripts.values())}")

    return list(transcripts.values())

"""
Collection of functions related to BED files

0-based
"""
# standard library imports
import collections
import sys
import time

# 3rd party library imports
# none

# local library imports
from g2gtools.exceptions import G2GBedError
from g2gtools.exceptions import G2GError
from g2gtools.exceptions import KeyboardInterruptError
from g2gtools.vci import VCIFile
import g2gtools.g2g_utils as g2g_utils

bed_fields = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'extra']
BEDRecord = collections.namedtuple('BEDRecord', bed_fields)

logger = g2g_utils.get_logger('g2gtools')


class BED(object):
    """
    Simple BED object for parsing and iterating BED files.

    Supports transparent gzip decompression.
    """

    def __init__(self, file_name: str):
        self.file_name: str = file_name
        self.current_line: str | None = None
        self.current_line_is_bed: bool = False
        self.current_record = None
        self.reader = g2g_utils.open_resource(file_name)
        self.n_items: int | None = None
        self.current_line_no: int = 0

    def __iter__(self):
        """
        Iteration.
        """
        return self

    def __next__(self):
        """
        Iteration.

        Raises:
            G2GBedError: When there is an improperly formatted BED file.
        """
        self.current_line = g2g_utils.s(self.reader.__next__())
        self.current_line_no += 1

        while self.current_line and len(self.current_line.strip()) == 0:
            self.current_line = self.reader.__next__()
            self.current_line_no += 1

        if self.current_line.startswith('track'):
            self.current_line = self.current_line.strip()
            self.current_line_is_bed = False
            self.current_record = None
            return None

        self.current_line_is_bed = True
        elem = self.current_line.strip().split('\t')

        if not self.n_items:
            self.n_items = len(elem)
        else:
            if self.n_items != len(elem):
                raise G2GBedError('Improperly formatted BED file')

        try:
            bed_data = {
                'chrom': elem[0],
                'start': int(elem[1]),
                'end': int(elem[2]),
                'name': elem[3] if self.n_items > 3 else None,
                'score': elem[4] if self.n_items > 4 else None,
                'strand': elem[5] if self.n_items > 5 else None,
                'extra': elem[6:] if self.n_items > 6 else None,
            }

            self.current_record = BEDRecord(**bed_data)
            return self.current_record
        except IndexError:
            raise G2GBedError(
                (
                    'Improperly formatted BED file, '
                    f'line number: {self.current_line_no}, '
                    f'line: {self.current_line}'
                )
            )
        except ValueError:
            raise G2GBedError(
                (
                    'Improperly formatted BED file, '
                    f'line number: {self.current_line_no}, '
                    f'line: {self.current_line}'
                )
            )


def convert_bed_file(
    vci_file: str | VCIFile,
    bed_file_name_in: str,
    bed_file_name_out: str = None,
    reverse: bool = False
) -> None:
    """
    Convert BED coordinates.

    Args
        vci_file: Name of the VCI file.
        bed_file_name_in: Input BED file to convert.
        bed_file_name_out: Name of output BED file.
        reverse: True to process VCI in reverse.
    """
    start = time.time()

    try:

        if not isinstance(vci_file, VCIFile):
            vci_file = g2g_utils.check_file(vci_file)
            vci_file = VCIFile(vci_file)
            logger.warning(f'Input VCI File: {vci_file.filename}')
            vci_file.parse(reverse)

        logger.warning(f'Input VCI File: {vci_file.filename}')

        bed_file_name_in = g2g_utils.check_file(bed_file_name_in)
        logger.warning(f'Input BED File: {bed_file_name_in}')

        if bed_file_name_out:
            bed_file_name_out = g2g_utils.check_file(bed_file_name_out, 'w')
            unmapped_file = f'{bed_file_name_out}.unmapped'
            bed_out = open(bed_file_name_out, 'w')
            bed_unmapped_file = open(unmapped_file, 'w')
            logger.info(f'Output BED File: {bed_file_name_out}')
            logger.info(f'Output UNMAPPED File: {unmapped_file}')
        else:
            input_dir, input_name = g2g_utils.get_dir_and_file(bed_file_name_in)
            unmapped_file = f'{input_name}.unmapped'
            unmapped_file = g2g_utils.check_file(unmapped_file, 'w')
            bed_out = sys.stdout
            bed_unmapped_file = open(unmapped_file, 'w')
            logger.info('Output BED File: stdout')
            logger.info(f'Output UNMAPPED File: {unmapped_file}')

        logger.info(f'VCI File is diploid: {vci_file.is_diploid()}')

        left_right = [''] if vci_file.is_haploid() else ['_L', '_R']

        logger.warning('Converting BED file...')

        bed_file = BED(bed_file_name_in)

        total = 0
        success = 0
        fail = 0

        for record in bed_file:
            logger.debug(f'ORIGINAL: {str(bed_file.current_line).strip()}')

            total += 1

            if total % 100000 == 0:
                logger.info(f'Processed {total:,} lines')

            for lr in left_right:
                seq_id = f'{record.chrom}{lr}'
                mappings = vci_file.find_mappings(
                    seq_id, record.start - 1, record.end
                )

                # unmapped
                if mappings is None:
                    logger.debug('\tFail due to no mappings')
                    bed_unmapped_file.write(bed_file.current_line)
                    fail += 0
                    continue
                else:
                    logger.debug(f'{len(mappings)} mappings found')

                success += 1
                start = mappings[0].to_start + 1
                end = mappings[-1].to_end
                elem = bed_file.current_line.rstrip().split('\t')

                logger.debug(f'({record.start-1}, {record.end})=>({start}, {end})')
                logger.debug(elem)

                elem[0] = seq_id
                elem[1] = start
                elem[2] = end

                temp_elem = '\t'.join(map(str, elem))
                logger.debug(f'     NEW: {temp_elem}')

                bed_out.write(f'{temp_elem}\n')

        if bed_file_name_out:
            bed_out.close()

        bed_unmapped_file.close()

        logger.warning(f'Converted {success:,} of {total:,} records')
        logger.warning('BED File converted')

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        raise G2GError(str(e))
    finally:
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warning(f'Time: {fmt_time}')

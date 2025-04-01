"""
A way to transform a Fasta file quickly.
"""
# standard library imports
from io import StringIO
import copy
import logging
import multiprocessing
import time

# 3rd party library imports
import pysam

# local library imports
from g2gtools.exceptions import G2GError
from g2gtools.exceptions import G2GFastaError
from g2gtools.exceptions import G2GRegionError
from g2gtools.exceptions import G2GValueError
from g2gtools.exceptions import KeyboardInterruptError
import g2gtools.fasta as fasta
import g2gtools.fasta_patch as fasta_patch
import g2gtools.g2g_utils as g2g_utils
import g2gtools.region as region
import g2gtools.vci as vci

logger = g2g_utils.get_logger('g2gtools')


class FastaTransformParams:
    """
    Parameters for transforming a FASTA file between genome assemblies.

    This class encapsulates all the parameters needed for the FASTA transformation
    process, including input/output file paths, regions, VCI (Variant Call
    Information) details, and other configuration options. It supports both
    coordinate transformation and optional sequence patching.

    Attributes:
        input_region (str | region.Region | None): Region to transform in the
            input file.
            For haploid genomes: "chr:start-end"
            For diploid genomes: "chr(_L|_R):start-end"
        input_file (str | None): Path to the input FASTA file.
        temp_dir (str | None): Directory for temporary files.
        full_file (bool): Whether to transform the entire file or just a region.
        output_file (str | None): Path to the output FASTA file.
        output_region (str | region.Region | None): Region identifier for the
            output file.
        output_header (str | fasta.FastaHeader | None): Header to use for the
             output FASTA sequence.
        vci_file (str | None): Path to the VCI file containing variants.
        vci_query (str | None): Query to filter variants from the VCI file.
        reverse (bool): Whether to reverse the transformation direction.
        gen_temp (bool): Whether to generate temporary files.
        offset (int): Offset to apply to coordinates.
        patch (bool): Whether to patch the sequence with variants.
        log_level (int | None): Logging level for the transformation process.
    """
    input_region: str | region.Region | None
    input_file: str | None
    temp_dir: str | None
    full_file: bool
    output_file: str | None
    output_region: str | region.Region | None
    output_header: str | fasta.FastaHeader | None
    vci_file: str | None
    vci_query: str | None
    reverse: bool
    gen_temp: bool
    offset: int
    patch: bool
    log_level: int | None

    def __init__(self) -> None:
        """
        Initialize a new FastaTransformParams object with default values.

        input region, must match the format of the input fasta file
        if fasta file is haploid, format is chr:start-end
        if fasta file is diploid, format is chr(_L|_R):start-end
        """
        self.input_region = None
        # original input fasta file
        self.input_file = None
        # temporary directory
        self.temp_dir = None
        self.full_file = False
        self.output_file = None
        self.output_region = None
        self.output_header = None
        # vci information
        self.vci_file = None
        self.vci_query = None
        self.reverse = False
        self.gen_temp = True
        # offset
        self.offset = 0
        # patch?
        self.patch = False
        self.log_level = logging.INFO

    def __str__(self) -> str:
        """
        Get a string representation of the transformation parameters.

        Returns:
            A string containing the input file, output file, location, offset,
            output region, output header, and VCI file.
        """
        return (
            f'Input: {self.input_file}\n\tOutput: {self.output_file}\n'
            f'\tLocation: {self.input_region}\n\tOffset: {self.offset}\n'
            f'\tOutput Region: {self.output_region}\n'
            f'\tOutput Header: {self.output_header}\n'
            f'\tVCI File: {self.vci_file}'
        )


class FastaTransformResult:
    """
    Result of a FASTA transformation operation.

    This class stores the results of transforming a FASTA file between genome
    assemblies, including the output file path and counts of different types
    of variants applied during the transformation.

    Attributes:
        output_file (str | None): Path to the output FASTA file.
        snp_count (int): Number of SNP (Single Nucleotide Polymorphism) variants applied.
        ins_count (int): Number of insertion variants applied.
        del_count (int): Number of deletion variants applied.
    """
    output_file: str | None
    snp_count: int
    ins_count: int
    del_count: int

    def __init__(self) -> None:
        """
        Initialize a new FastaTransformResult object with default values.
        """
        self.output_file = None
        self.snp_count = 0
        self.ins_count = 0
        self.del_count = 0

    def __str__(self) -> str:
        """
        Get a string representation of the transformation result.

        Returns:
            A string containing the output file path.
        """
        return f'File: {self.output_file}'


def process_piece(
    params: FastaTransformParams,
) -> FastaTransformResult:
    """
    Process the 'piece' of information.

    Args:
        params: The parameters dictating what piece to process.

    Returns:
        The results of the processing.
    """
    logger = g2g_utils.configure_logging('g2gtools', params.log_level)
    # TODO: check patch and make sure that we are using the same type of info
    logger.debug(f'params:input_region={params.input_region}')
    logger.debug(f'params:input_file={params.input_file}')
    logger.debug(f'params:output_file={params.output_file}')
    logger.debug(f'params:output_region={params.output_region}')
    logger.debug(f'params:output_header={params.output_header}')
    logger.debug(f'params:vci_query={params.vci_query}')

    fasta_transform_result = FastaTransformResult()
    fasta_transform_result.output_file = params.output_file

    try:
        fasta_file = fasta.FastaFile(params.input_file)
        prefix = g2g_utils.location_to_filestring(params.input_region)
        tmp_fasta = g2g_utils.gen_file_name(
            prefix=f'{prefix}_',
            extension='fa',
            append_time=False,
            output_dir=params.temp_dir,
        )
        logger.debug(f'tmp_fasta={tmp_fasta}')

        working = open(tmp_fasta, 'w')
        working.write(f'>{params.input_region.seq_id}\n')
        logger.debug(f'Fasta Fetch {params.input_region}')

        sequence = fasta_file.fetch(
            params.input_region.seq_id,
            params.input_region.start - 1,
            params.input_region.end,
        )

        if len(sequence) < 50:
            logger.debug(f'Fasta Fetch = {sequence}')
        else:
            logger.debug(f'Fasta Fetch = {sequence[:25]} ... {sequence[-25:]}')

        g2g_utils.write_sequence(sequence, working)
        working.close()

        if params.patch:
            logger.info('####################################################')
            logger.info('################## PATCHING ########################')
            logger.info('####################################################')

            p = fasta_patch.FastaPatchParams()
            p.input_region = params.input_region

            # original input fasta file
            p.input_file = tmp_fasta

            # temporary directory
            p.temp_dir = params.temp_dir

            p.output_file = tmp_fasta
            p.output_region = params.input_region
            f_desc = (
                f'{params.input_region.seq_id} '
                f'{params.input_region.seq_id}:'
                f'{params.input_region.start}-'
                f'{params.input_region.end}'
            )
            p.output_header = fasta.FastaHeader(
                str(params.input_region.seq_id), f_desc
            )

            # vci information
            p.vci_file = params.vci_file
            p.vci_query = params.vci_query
            p.reverse = params.reverse

            p.gen_temp = False

            # offset
            p.offset = params.offset

            p1 = fasta_patch.process_piece(p)
            logger.debug(f'\npatch output = {p1.output_file}\n')
            fasta_file = fasta.FastaFile(p1.output_file)
            fasta_transform_result.snp_count = p1.count
        else:
            fasta_file = fasta.FastaFile(tmp_fasta)

        offset = params.offset
        reverse = params.reverse

        logger.info(f'Transforming {params.input_region}...')
        reg = region.parse_region(params.vci_query)
        logger.info(f'Finding VCI mappings for {reg}')

        logger.info('PARSING VCI FILE...')
        vci_file = vci.VCIFile(params.vci_file, seq_ids=[reg.seq_id])
        vci_file.parse(reverse=reverse)
        logger.info('DONE PARSING!')

        mappings = vci_file.find_mappings(
            reg.seq_id,
            params.output_region.start - 1,
            params.output_region.end,
        )
        fasta_out = open(params.output_file, 'w')

        if mappings is None:
            logger.info('This region was deleted')
            logger.info('TODO: dump the fasta sequence here')
            logger.info(
                'output_region.seq_id=' f'{params.output_region.seq_id}'
            )
            logger.info(f'output_file={params.output_file}')

            if params.full_file:
                out_header = (
                    f'>{params.output_region.seq_id} '
                    f'{params.input_region.seq_id}:'
                    f'{params.input_region.start}-'
                    f'{params.input_region.end} from|'
                    f'{params.input_region.seq_id}:'
                    f'{params.input_region.start}-'
                    f'{params.input_region.end}\n'
                )
            else:
                out_header = (
                    f'>{params.input_region.seq_id}:'
                    f'{params.input_region.start}-'
                    f'{params.input_region.end} from|'
                    f'{params.input_region.seq_id}:'
                    f'{params.input_region.start}-'
                    f'{params.input_region.end}\n'
                )

            fasta_out.write(out_header)

            partial_seq = fasta_file.fetch(
                params.input_region.seq_id,
                params.input_region.start - 1,
                params.input_region.end,
            )

            if len(partial_seq) < 50:
                logger.debug(f'Fasta Fetch = {partial_seq}')
            else:
                logger.debug(
                    f'Fasta Fetch = {partial_seq[:25]}...{partial_seq[-25:]}'
                )

            g2g_utils.write_sequence(partial_seq, fasta_out)

            return fasta_transform_result

        first = True
        fasta_transform_result.ins_count = 0
        fasta_transform_result.del_count = 0
        new_start_pos = mappings[0].to_start
        # new_end_pos = mappings[-1].to_end
        last_pos = 0
        new_sequence = StringIO()
        start = mappings[0].from_start
        found = False
        new_sequence_len = 0

        logger.debug(f'new_start_pos={new_start_pos}')
        logger.debug(f'Setting start to {start} (mappings[0].from_start)')
        logger.debug(f'start={start}')
        logger.debug(f'last_pos={last_pos}')
        logger.debug(f'VCI Fetch {params.vci_query}')

        for line in vci_file.fetch(params.vci_query, parser=pysam.asTuple()):
            # CHROM POS ANCHOR DEL INS FRAG
            if line[5] == '.':
                continue

            logger.debug(f'LINE: {line}')

            found = True
            shared_bases = line[2]
            deleted_bases = line[3 if not reverse else 4]

            if deleted_bases != '.':
                deleted_bases_length = len(deleted_bases)
            else:
                deleted_bases_length = 0

            inserted_bases = line[4 if not reverse else 3]

            if inserted_bases != '.':
                inserted_bases_length = len(inserted_bases)
            else:
                inserted_bases_length = 0

            fragment_size = int(line[5])

            if first:
                # fragment_size = (current_pos + shared_bases_length) -
                #                 (previous_pos + previous_shared_length +
                #                               + previous_inserted_length)
                logger.debug('First result in query...')
                logger.debug(f'Adjusting last_pos from {last_pos} to {start}')
                last_pos = start

                old_fragment_size = fragment_size
                fragment_size = (int(line[1]) + len(shared_bases)) - (
                    last_pos + 1
                )
                logger.debug(
                    f'Adjusting fragment_size from {old_fragment_size} '
                    f'to {fragment_size}'
                )

                first = False

                if fragment_size < 0:
                    continue

            logger.debug(f'last_pos={last_pos}')
            logger.debug(f'offset={offset}')
            logger.debug(f'fragment_size={fragment_size}')
            logger.debug(
                f'extracting... {last_pos - offset}-'
                f'{last_pos + fragment_size - offset}'
            )
            partial_seq = fasta_file.fetch(
                params.input_region.seq_id,
                last_pos - offset,
                last_pos + fragment_size - offset,
            )

            if len(partial_seq) < 50:
                logger.debug(f'Fasta Fetch = {partial_seq}')
            else:
                logger.debug(
                    f'Fasta Fetch = {partial_seq[:25]} ... {partial_seq[-25:]}'
                )

            new_sequence.write(partial_seq)
            new_sequence_len += fragment_size

            if inserted_bases_length > 0:
                # insertion
                logger.debug('INSERTION')
                logger.debug(f'Adding {inserted_bases}')

                new_sequence.write(inserted_bases)
                fasta_transform_result.ins_count += inserted_bases_length
                new_sequence_len += inserted_bases_length

            if deleted_bases_length > 0:
                # deletion
                logger.debug('DELETION')

                last_pos += deleted_bases_length
                fasta_transform_result.del_count += deleted_bases_length

            last_pos += fragment_size

        if found:
            logger.debug('Fetching last bit of sequence')

            if last_pos >= params.input_region.end:
                logger.debug('Nothing to fetch... done')
            else:
                logger.debug(
                    f'Fasta Fetch {params.input_region.seq_id}:'
                    f'{last_pos - offset}-'
                    f'{params.input_region.end - offset} '
                    '(0-based)'
                )

                partial_seq = fasta_file.fetch(
                    params.input_region.seq_id,
                    last_pos - offset,
                    params.input_region.end - offset,
                )

                if len(partial_seq) < 50:
                    logger.debug(f'Fasta Fetch = {partial_seq}')
                else:
                    logger.debug(
                        f'Fasta Fetch = {partial_seq[:25]} ... {partial_seq[-25:]}'
                    )

                new_sequence.write(partial_seq)
                new_sequence_len += params.input_region.end - last_pos
        else:
            logger.debug('NOTHING FOUND IN REGION')
            logger.debug('Fetching ONLY bit of sequence')
            logger.debug(f'start={start}')
            logger.debug(f'offset={offset}')
            logger.debug(
                'params.input_region.end=' f'{params.input_region.end}'
            )

            logger.debug(
                f'Fasta Fetch {params.input_region.seq_id}:'
                f'{start - offset}-'
                f'{params.input_region.end - offset} (0-based)'
            )
            partial_seq = fasta_file.fetch(
                params.input_region.seq_id,
                start - offset,
                params.input_region.end - offset,
            )
            new_sequence.write(partial_seq)
            new_sequence_len += len(partial_seq)

        logger.debug(f'the new_sequence_len={new_sequence_len}')

        if params.full_file:
            logger.debug('FULL FILE')
            out_header = (
                f'>{params.output_region.seq_id} '
                f'{params.input_region.seq_id}:'
                f'{new_start_pos + 1}-{new_start_pos + new_sequence_len} '
                f'from|{params.input_region.seq_id}:'
                f'{params.input_region.start}-'
                f'{params.input_region.end}\n'
            )
        else:
            logger.debug('NOT FULL FILE')
            out_header = (
                f'>{params.input_region.seq_id}:'
                f'{new_start_pos + 1}-{new_start_pos + new_sequence_len} '
                f'from|{params.input_region.seq_id}:'
                f'{params.input_region.start}-'
                f'{params.input_region.end}\n'
            )

        logger.debug(f'WRITING HEADER: {out_header}')
        fasta_out.write(out_header)
        g2g_utils.write_sequence(new_sequence.getvalue(), fasta_out)
        fasta_out.close()
    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except G2GRegionError as le:
        logger.debug(f'Unable to parse (1) location, {le}')
        raise le
    except G2GValueError as e:
        logger.debug(f'Unable to parse (2) location, {e}')
        raise e
    except G2GFastaError as e:
        logger.debug(f'Unable to parse (3) location, {e}')
        raise e
    except TypeError as e:
        logger.debug(f'Unable to parse (4) location, {e}')
        raise e
    except Exception as e:
        logger.debug(f'Unable to parse (5) location, {e}')
        raise e

    logger.info(f'Transforming complete: {params.input_region}')

    g2g_utils.delete_index_files(tmp_fasta)
    g2g_utils.delete_file(tmp_fasta)

    return fasta_transform_result


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    Args:
        args: The arguments to process_piece.

    Returns:
        The same as process_piece
    """
    return process_piece(*args)


def prepare_fasta_transform(filename_output: str) -> str:
    """
    Get the complete file path and make sure it is ok to write.

    Args:
        filename_output: The file name.

    Returns:
        The complete path to the file.

    Raises:
        G2GValueError: If the fasta name doesn't end correctly.
    """
    filename_output = g2g_utils.check_file(filename_output, 'w')

    new_filename_output = filename_output

    # let's figure out what our output name will be
    if new_filename_output.lower().endswith('.gz'):
        # strip off .gz
        new_filename_output = new_filename_output[:-3]

    if not filename_output.lower().endswith('.fa'):
        raise G2GValueError(
            'Expecting output filename extension to be either ".fa.gz" or ".fa"'
        )

    g2g_utils.delete_index_files(new_filename_output)

    return new_filename_output


def process(
    filename_fasta: str,
    filename_vci: str,
    regions: list[region.Region] | None = None,
    filename_output: str | None = None,
    bgzip: bool | None = False,
    reverse: bool | None = False,
    num_processes: int | None = None,
    also_patch: bool | None = False,
) -> None:
    """
    Patch a Fasta file by replacing the bases where the SNPs are located in
    the VCI file.

    Args
        filename_fasta: Input Fasta file to convert.
        filename_vci: Name of the VCI file or a VCIFile object.
        regions: A list of Regions to process.
        filename_output: Name of output Fasta file.
        bgzip: True to compress and index file.
        reverse: Reverse the insertions and deletions in VCI file.
        num_processes: The number of processes to use.
        also_patch: True to also patch the Fasta file.
    """
    start = time.time()
    dump_fasta = False
    temp_directory = g2g_utils.create_temp_dir('transform_', dir='.')
    filename_fasta = g2g_utils.check_file(filename_fasta)
    filename_vci = g2g_utils.check_file(filename_vci)

    if not num_processes:
        num_processes = multiprocessing.cpu_count()
    else:
        if num_processes <= 0:
            num_processes = 1

    logger.warning(f'Input VCI File: {filename_vci}')
    logger.warning(f'Input Fasta File: {filename_fasta}')
    logger.debug(f'Number of processes: {num_processes}')
    logger.debug(f'Temp directory: {temp_directory}')

    try:
        if filename_output:
            filename_output = g2g_utils.check_file(filename_output, 'w')

            if not regions:
                filename_output = prepare_fasta_transform(filename_output)
                logger.warning(f'Output Fasta File: {filename_output}')
            else:
                if bgzip:
                    if filename_output.lower().endswith(('.fa', '.fasta')):
                        logger.warning(
                            f'Output Fasta File: {filename_output}.gz'
                        )
                    elif filename_output.lower().endswith('.gz'):
                        logger.warning(f'Output Fasta File: {filename_output}')
                        filename_output = filename_output[:-3]
                else:
                    logger.warning(f'Output Fasta File: {filename_output}')
        else:
            filename_output = g2g_utils.gen_file_name(
                extension='fa', append_time=False, output_dir=temp_directory
            )
            dump_fasta = True
            logger.debug(f'Temporary Fasta File: {filename_output}')

        fasta_file = fasta.FastaFile(filename_fasta)
        vci_file = vci.VCIFile(filename_vci)

        if fasta_file.is_diploid() and vci_file.is_haploid():
            raise G2GFastaError(
                'Haploid VCI file and diploid Fasta file combination '
                'is not currently supported for transform'
            )

        full_file = True

        if regions:
            full_file = False
            if len(regions) > 5:
                regions_str = ', '.join(reg for reg in map(str, regions[:5]))
                logger.warning(f'Regions: {regions_str} (showing 1st five)')
            else:
                regions_str = ', '.join(reg for reg in map(str, regions))
                logger.warning(f'Regions: {regions_str}')

        else:
            logger.debug('No regions, calculating from Fasta file...')
            regions = []
            for chrom in fasta_file.references:
                regions.append(
                    region.Region(
                        chrom, 1, fasta_file.get_reference_length(chrom)
                    )
                )

        all_params = []
        for reg in regions:
            logger.debug(f'reg={reg}')
            logger.debug(f'reg.original_base={reg.original_base}')

            params = FastaTransformParams()
            params.input_region = reg
            params.input_file = filename_fasta
            params.temp_dir = temp_directory
            params.vci_file = filename_vci
            params.reverse = reverse
            params.offset = 0 if reg.start <= 1 else reg.start - 1
            params.patch = also_patch
            params.full_file = full_file
            params.log_level = logger.level

            if fasta_file.is_haploid() and vci_file.is_diploid():
                # logger.info('*** Experimental ***')
                # logger.info('*** HAPLOID FASTA and DIPLOID VCI ***')
                logger.info('Fasta file is haploid and VCI file is diploid')
                prefix = g2g_utils.location_to_filestring(reg)
                params.output_file = g2g_utils.gen_file_name(
                    prefix=f'{prefix}_L',
                    extension='fa',
                    output_dir=temp_directory,
                    append_time=False,
                )

                if full_file:
                    params.output_region = region.Region(
                        f'{reg.seq_id}_L', reg.start, reg.end
                    )
                    params.output_header = fasta.FastaHeader(
                        f'{reg.seq_id}_L',
                        f'{reg.seq_id}:{reg.start}-{reg.end}',
                    )
                else:
                    params.output_region = region.Region(
                        f'{reg.seq_id}_L', reg.start, reg.end
                    )
                    params.output_header = fasta.FastaHeader(
                        f'{reg.seq_id}_L',
                        f'{reg.seq_id}_L:{reg.start}-{reg.end}',
                    )

                params.vci_query = f'{reg.seq_id}_L:{reg.start}-{reg.end}'
                all_params.append(params)

                params_r = copy.deepcopy(params)
                prefix_r = g2g_utils.location_to_filestring(reg)
                params_r.output_file = g2g_utils.gen_file_name(
                    prefix=f'{prefix_r}_R',
                    extension='fa',
                    output_dir=temp_directory,
                    append_time=False,
                )

                if full_file:
                    params_r.output_region = region.Region(
                        f'{reg.seq_id}_R', reg.start, reg.end
                    )
                    params_r.output_header = fasta.FastaHeader(
                        f'{reg.seq_id}_R',
                        f'{reg.seq_id}:{reg.start}-{reg.end}',
                    )
                else:
                    params_r.output_region = region.Region(
                        f'{reg.seq_id}_R', reg.start, reg.end
                    )
                    params_r.output_header = fasta.FastaHeader(
                        f'{reg.seq_id}_R',
                        f'{reg.seq_id}_R:{reg.start}-{reg.end}',
                    )

                params_r.vci_query = f'{reg.seq_id}_R:{reg.start}-{reg.end}'
                all_params.append(params_r)
            else:
                fmt_string = 'HAPLOID' if vci_file.is_haploid() else 'DIPLOID'
                logger.debug(f'VCI file and Fasta file are both: {fmt_string}')
                params.output_file = g2g_utils.gen_file_name(
                    prefix=g2g_utils.location_to_filestring(reg),
                    extension='fa',
                    output_dir=temp_directory,
                    append_time=False,
                )

                if full_file:
                    params.output_region = region.Region(
                        reg.seq_id, reg.start, reg.end
                    )
                    params.output_header = fasta.FastaHeader(
                        reg.seq_id,
                        f'{reg.seq_id}:{reg.start}-{reg.end}',
                    )
                else:
                    params.output_region = region.Region(
                        reg.seq_id, reg.start, reg.end
                    )
                    params.output_header = fasta.FastaHeader(
                        f'{reg.seq_id} ' '{reg.seq_id}:{reg.start}-{reg.end}',
                        None,
                    )

                params.vci_query = f'{reg.seq_id}:{reg.start}-{reg.end}'
                all_params.append(params)

        pool = multiprocessing.Pool(num_processes)
        args = zip(all_params)
        results = pool.map(wrapper, args)

        # parse results
        total_snp = 0
        total_ins = 0
        total_del = 0

        files = []

        for c in results:
            if c is not None:
                files.append(c.output_file)
                total_snp += c.snp_count
                total_ins += c.ins_count
                total_del += c.del_count

        if also_patch:
            logger.warning(f'Processed {total_snp:,} SNPs')

        logger.warning(f'Processed {total_ins:,} insertions')
        logger.warning(f'Processed {total_del:,} deletions')

        logger.debug('Temp files created, copy to master temp file and delete')
        g2g_utils.concatenate_files(files, filename_output, False)

        if dump_fasta:
            g2g_utils.dump_file_contents(filename_output)
        else:
            # move temp to final destination
            if bgzip:
                # remove the fai
                logger.debug(f'Removing the FAI index for {filename_output}')
                g2g_utils.delete_index_files(filename_output)

                logger.warning('Compressing and indexing...')

                if filename_output.lower().endswith(('.fa', '.fasta')):
                    g2g_utils.bgzip_and_index_file(
                        filename_output,
                        f'{filename_output}.gz',
                        delete_original=True,
                        file_format='fa',
                    )
                elif filename_output.lower().endswith('.gz'):
                    g2g_utils.bgzip_and_index_file(
                        filename_output,
                        filename_output,
                        delete_original=True,
                        file_format='fa',
                    )

        logger.warning('Fasta File transformed')

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        raise G2GError(str(e))
    finally:
        # clean up the temporary files
        g2g_utils.delete_dir(temp_directory)
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warning(f'Time: {fmt_time}')

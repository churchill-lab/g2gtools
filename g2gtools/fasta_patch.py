#
# A way to transform a Fasta file quickly
#

# standard library imports
import copy
import mmap
import multiprocessing
import time

# 3rd party library imports
import pysam

# local library imports
from .exceptions import G2GError
from .exceptions import G2GFastaError
from .exceptions import G2GValueError
from .exceptions import KeyboardInterruptError
from . import fasta
from . import g2g
from . import g2g_utils
from . import vci


class FastaPatchParams(object):
    def __init__(self):
        # input region, must match the format of the input fasta file
        # if fasta file is haploid, format is chr:start-end
        # if fasta file is diploid, format is chr(_L|_R):start-end
        self.input_region = None

        # original input fasta file
        self.input_file = None

        # temporary directory
        self.temp_dir = None

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

        self.debug_level = 0

    def __str__(self):
        return (
            f"Input: {self.input_file}\nOutput: {self.output_file}\n"
            f"Location: {self.input_region}\nOffset: {self.offset}"
        )


class FastaPatchResult(object):
    def __init__(self):
        self.output_file = None
        self.count = 0

    def __str__(self):
        return f"File: {self.output_file}"


def process_piece(fasta_patch_params: FastaPatchParams) -> FastaPatchResult:
    """
    Process the "piece" of information.

    Args:
        fasta_patch_params: The parameters dictating what piece to process.

    Returns:
        The results of the processing.
    """
    logger = g2g.get_logger(fasta_patch_params.debug_level)

    logger.debug(f"params:input_region = {fasta_patch_params.input_region}")
    logger.debug(f"params:input_file = {fasta_patch_params.input_file}")
    logger.debug(f"params:temp_dir = {fasta_patch_params.temp_dir}")
    logger.debug(f"params:output_file = {fasta_patch_params.output_file}")
    logger.debug(f"params:output_region = {fasta_patch_params.output_region}")
    logger.debug(f"params:output_header = {fasta_patch_params.output_header}")
    logger.debug(f"params:vci_file = {fasta_patch_params.vci_file}")
    logger.debug(f"params:vci_query = {fasta_patch_params.vci_query}")
    logger.debug(f"params:reverse = {fasta_patch_params.reverse}")
    logger.debug(f"params:gen_temp = {fasta_patch_params.gen_temp}")
    logger.debug(f"params:offset = {fasta_patch_params.offset}")

    fasta_patch_result = FastaPatchResult()
    fasta_patch_result.output_file = fasta_patch_params.output_file

    byte_start, byte_end, byte_len_seq = 0, 0, 0
    line = None
    mm = None

    try:
        fasta_file = fasta.FastaFile(fasta_patch_params.input_file)
        vci_file = vci.VCIFile(fasta_patch_params.vci_file)

        if fasta_patch_params.gen_temp:
            prefix = g2g_utils.location_to_filestring(
                fasta_patch_params.input_region
            )
            tmp_fasta = g2g_utils.gen_file_name(
                prefix=f"{prefix}_", extension="fa", append_time=False,
                output_dir=fasta_patch_params.temp_dir
            )
            logger.debug(f"tmp_fasta={tmp_fasta}")

            working = open(tmp_fasta, "w")
            if fasta_patch_params.output_header.description:
                fasta_header_string = (
                    f">{fasta_patch_params.output_header.id} "
                    f"{fasta_patch_params.output_header.description}"
                )
            else:
                fasta_header_string = f">{fasta_patch_params.output_header.id}"
            working.write(f"{fasta_header_string}\n")

            logger.debug(f"Fasta Fetch {fasta_patch_params.input_region}")
            logger.debug(f"Fasta Fetch {fasta_patch_params.input_region.start}")
            sequence = fasta_file.fetch(
                fasta_patch_params.input_region.seq_id,
                fasta_patch_params.input_region.start - 1,
                fasta_patch_params.input_region.end
            )

            if len(sequence) < 50:
                logger.debug(f"Fasta Fetch = {sequence}")
            else:
                logger.debug(
                    f"Fasta Fetch = {sequence[:25]} ... {sequence[-25:]}"
                )
            g2g_utils.write_sequence(sequence, working)
            working.close()

            fasta.FastaFile(tmp_fasta)
            fasta_patch_result.output_file = tmp_fasta
            fasta_index = fasta.FAI(tmp_fasta)
        else:
            fasta_index = fasta.FAI(fasta_patch_params.input_file)
            tmp_fasta = fasta_patch_params.input_file

        offset = fasta_patch_params.offset
        reverse = fasta_patch_params.reverse

        logger.info(f"Processing {fasta_patch_params.output_header.id}...")

        fd_l = open(tmp_fasta, "r+b")
        mm = mmap.mmap(fd_l.fileno(), 0)

        logger.debug(f"VCI Fetch {fasta_patch_params.vci_query}")
        for line in vci_file.fetch(
                fasta_patch_params.vci_query, parser=pysam.asTuple()
        ):
            # CHROM POS ANCHOR DEL INS FRAG
            if line[5] != '.':
                continue

            position = int(line[1]) - offset
            deleted_bases = line[3 if not reverse else 4]
            inserted_bases = line[4 if not reverse else 3]

            logger.debug(f"int(line[1])={int(line[1])}")
            logger.debug(f"position={position}")
            logger.debug(f"offset={offset}")

            fasta_patch_result.count += 1

            byte_start, byte_end, byte_len_seq = fasta_index.get_pos(
                fasta_patch_params.output_region.seq_id, position - 1, position
            )

            logger.debug("LINE: {line}")
            logger.debug("WAS {chr(mm[byte_start])}")
            logger.debug(
                f"Patching {fasta_patch_params.output_region.seq_id}:"
                f"{position - 1} from {deleted_bases} to {inserted_bases}"
            )
            mm[byte_start] = ord(inserted_bases)
    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except ValueError:
        logger.debug(f"No SNPS found in region {fasta_patch_params.vci_query}")
    except Exception as e:
        logger.debug(e)
        logger.debug(f"{byte_start}, {byte_end}, {byte_len_seq}")
        logger.debug(line)

    if mm:
        mm.flush()
        mm.close()

    return fasta_patch_result


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    Args:
        args: The arguments to process_piece.

    Returns:
        The same as process_piece
    """
    return process_piece(*args)


def prepare_fasta_patch(filename_output: str) -> str:
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
            "Expecting output filename extension to be '.fa.gz' or '.fa'"
        )

    g2g_utils.delete_index_files(new_filename_output)

    return new_filename_output


def process(
        filename_fasta: str,
        filename_vci: str,
        regions: list[g2g.Region] | None = None,
        filename_output: str | None = None,
        bgzip: bool = False,
        reverse: bool = False,
        num_processes: int | None = None,
        debug_level: int = 0
) -> None:
    """
    Patch a Fasta file by replacing the bases where the SNPs are located as
    specified in the VCI file created from the VCF file.

    Args:
        filename_fasta: The name of the Fasta file.
        filename_vci: The name of the VCI file.
        regions: A specific list of Regions to use or None for the entire file.
        filename_output: The of the output file.
        bgzip: True to compress and index file.
        reverse: Reverse the insertions and deletions in VCI file.
        num_processes: The number of processes to use.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    start = time.time()
    dump_fasta = False
    temp_directory = g2g_utils.create_temp_dir('patch_', dir='.')
    filename_fasta = g2g_utils.check_file(filename_fasta)
    filename_vci = g2g_utils.check_file(filename_vci)

    if not num_processes:
        num_processes = multiprocessing.cpu_count()
    else:
        if num_processes <= 0:
            num_processes = 1

    logger = g2g.get_logger(debug_level)
    logger.warn(f"Processes: {num_processes}")
    logger.warn(f"Input VCI File: {filename_vci}")
    logger.warn(f"Input Fasta File: {filename_fasta}")
    logger.debug(f"Temp directory: {temp_directory}")

    try:
        if filename_output:
            filename_output = g2g_utils.check_file(filename_output, 'w')

            if not regions:
                filename_output = prepare_fasta_patch(filename_output)
                logger.warn(f"Output Fasta File: {filename_output}")
            else:
                if bgzip:
                    if filename_output.lower().endswith((".fa", ".fasta")):
                        logger.warn(f"Output Fasta File: {filename_output}.gz")
                    elif filename_output.lower().endswith(".gz"):
                        logger.warn(f"Output Fasta File: {filename_output}")
                        filename_output = filename_output[:-3]
                else:
                    logger.warn(f"Output Fasta File: {filename_output}")
        else:
            filename_output = g2g_utils.gen_file_name(
                extension="fa", append_time=False, output_dir=temp_directory
            )
            dump_fasta = True
            logger.debug(f"Temporary fasta file: {filename_output}")

        fasta_file = fasta.FastaFile(filename_fasta)
        vci_file = vci.VCIFile(filename_vci)

        if fasta_file.is_diploid():
            raise G2GFastaError(
                "Diploid Fasta files are not currently supported for patch"
            )

        full_file = True

        if regions:
            full_file = False
            if len(regions) > 5:
                regions_str = ", ".join(reg for reg in map(str, regions[:5]))
                logger.warn(f"Regions: {regions_str} (showing 1st five)")
            else:
                regions_str = ", ".join(reg for reg in map(str, regions))
                logger.warn(f"Regions: {regions_str}")

        else:
            logger.debug("No regions, calculating from Fasta file...")
            regions = []
            for chrom in fasta_file.references:
                regions.append(
                    g2g.Region(chrom, 1, fasta_file.get_reference_length(chrom))
                )

        all_params = []
        for region in regions:
            logger.debug(f"region={region}")
            logger.debug(f"region.original_base={region.original_base}")

            params = FastaPatchParams()
            params.input_region = region
            params.input_file = filename_fasta
            params.temp_dir = temp_directory
            params.vci_file = filename_vci
            params.reverse = reverse
            params.offset = 0 if region.start <= 1 else region.start - 1
            params.debug_level = debug_level

            if vci_file.is_haploid():
                logger.debug("VCI File is haploid")
                params.output_file = g2g_utils.gen_file_name(
                    prefix=g2g_utils.location_to_filestring(region),
                    extension="fa", output_dir=temp_directory,
                    append_time=False
                )

                if full_file:
                    params.output_region = g2g.Region(
                        region.seq_id, region.start, region.end
                    )
                    params.output_header = fasta.FastaHeader(
                        region.seq_id,
                        f"{region.seq_id}:{region.start}-{region.end}"
                    )
                else:
                    params.output_region = g2g.Region(
                        f"{region.seq_id}:{region.start}-{region.end}",
                        region.start, region.end
                    )
                    params.output_header = fasta.FastaHeader(
                        f"{region.seq_id}:{region.start}-{region.end}",
                        None
                    )

                params.vci_query = (
                    f"{region.seq_id}:{region.start}-{region.end}"
                )
                all_params.append(params)
            else:
                logger.debug("VCI File is diploid")
                params.output_file = g2g_utils.gen_file_name(
                    prefix=g2g_utils.location_to_filestring(region)+"_L",
                    extension="fa", output_dir=temp_directory,
                    append_time=False
                )

                if full_file:
                    params.output_region = g2g.Region(
                        f"{region.seq_id}_L", region.start, region.end
                    )
                    params.output_header = fasta.FastaHeader(
                        f"{region.seq_id}_L",
                        f"{region.seq_id}:{region.start}-{region.end}"
                    )
                else:
                    params.output_region = g2g.Region(
                        f"{region.seq_id}_L:{region.start}-{region.end}",
                        region.start, region.end
                    )
                    params.output_header = fasta.FastaHeader(
                        f"{region.seq_id}_L:{region.start}-{region.end}",
                        None
                    )

                params.vci_query = (
                    f"{region.seq_id}_L:{region.start}-{region.end}"
                )
                all_params.append(params)

                params_r = copy.deepcopy(params)
                params_r.output_file = g2g_utils.gen_file_name(
                    prefix=g2g_utils.location_to_filestring(region)+"_R",
                    extension="fa", output_dir=temp_directory, append_time=False
                )

                if full_file:
                    params_r.output_region = g2g.Region(
                        f"{region.seq_id}_R", region.start, region.end
                    )
                    params_r.output_header = fasta.FastaHeader(
                        region.seq_id+"_R",
                        f"{region.seq_id}:{region.start}-{region.end}"
                    )
                else:
                    params_r.output_region = g2g.Region(
                        f"{region.seq_id}_R:{region.start}-{region.end}",
                        region.start, region.end
                    )
                    params_r.output_header = fasta.FastaHeader(
                        f"{region.seq_id}_R:{region.start}-{region.end}",
                        None
                    )

                params_r.vci_query = (
                    f"{region.seq_id}_R:{region.start}-{region.end}"
                )
                all_params.append(params_r)

        args = zip(all_params)
        pool = multiprocessing.Pool(num_processes)
        # changed from imap_unordered to imap to keep original seq_id order
        results = pool.imap(wrapper, args)

        # parse results
        total = 0
        mode = 'wb'

        for c in results:
            if c is not None:
                total += c.count
                g2g_utils.concatenate_files(
                    [c.output_file], filename_output, False, mode
                )
                g2g_utils.delete_file(c.output_file)
                g2g_utils.delete_index_files(c.output_file)
                mode = "ab"

        logger.warn(f"Patched {total:,} SNPs total")

        if dump_fasta:
            g2g_utils.dump_file_contents(filename_output)
        else:
            # move temp to final destination
            if bgzip:
                # remove the fai
                logger.debug(f"Removing the FAI index for {filename_output}")
                g2g_utils.delete_index_files(filename_output)

                logger.warn("Compressing and indexing...")

                if filename_output.lower().endswith((".fa", ".fasta")):
                    g2g_utils.bgzip_and_index_file(
                        filename_output, f"{filename_output}.gz",
                        delete_original=True, file_format="fa"
                    )
                elif filename_output.lower().endswith(".gz"):
                    g2g_utils.bgzip_and_index_file(
                        filename_output, filename_output,
                        delete_original=True, file_format="fa"
                    )

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        raise G2GError(str(e))
    finally:
        # clean up the temporary files
        g2g_utils.delete_dir(temp_directory)
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warn(f"Patch complete: {fmt_time}")

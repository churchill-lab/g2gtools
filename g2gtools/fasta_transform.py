#
# A way to transform a Fasta file quickly
#

# standard library imports
from io import StringIO
import copy
import multiprocessing
import os
import time

# 3rd party library imports
import pysam

# local library imports
from .exceptions import G2GError
from .exceptions import G2GFastaError
from .exceptions import G2GRegionError
from .exceptions import G2GValueError
from .exceptions import KeyboardInterruptError
from . import fasta
from . import fasta_patch
from . import g2g
from . import g2g_utils
from . import vci


class FastaTransformParams(object):
    def __init__(self):
        # input region, must match the format of the input fasta file
        # if fasta file is haploid, format is chr:start-end
        # if fasta file is diploid, format is chr(_L|_R):start-end
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

        self.debug_level = 0

    def __str__(self):
        return (
            f"Input: {self.input_file}\n\tOutput: {self.output_file}\n"
            f"\tLocation: {self.input_region}\n\tOffset: {self.offset}\n"
            f"\tOutput Region: {self.output_region}\n"
            f"\tOutput Header: {self.output_header}\n"
            f"\tVCI File: {self.vci_file}"
        )


class FastaTransformResult(object):
    def __init__(self):
        self.output_file = None
        self.snp_count = 0
        self.ins_count = 0
        self.del_count = 0

    def __str__(self):
        return f"File: {self.output_file}"


def process_piece(fasta_transform_params):
    """
    """
    # todo: check with patch and make sure that we are using the same type of info
    start_time = time.time()
    logger = g2g.get_logger(fasta_transform_params.debug_level)
    fasta_transform_result = FastaTransformResult()
    fasta_transform_result.output_file = fasta_transform_params.output_file

    try:
        fasta_file = fasta.FastaFile(fasta_transform_params.input_file)

        logger.debug(f"fasta_transform_params.input_region={fasta_transform_params.input_region}")
        logger.debug(f"fasta_transform_params.input_file={fasta_transform_params.input_file}")
        logger.debug(f"fasta_transform_params.output_file={fasta_transform_params.output_file}")
        logger.debug(f"fasta_transform_params.output_region={fasta_transform_params.output_region}")
        logger.debug(f"fasta_transform_params.output_header={fasta_transform_params.output_header}")
        logger.debug(f"fasta_transform_params.vci_query={fasta_transform_params.vci_query}")

        tmp_fasta = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(fasta_transform_params.input_region) + "_", extension="fa", append_time=False, output_dir=fasta_transform_params.temp_dir)
        logger.debug(f"tmp_fasta={tmp_fasta}")

        working = open(tmp_fasta, "w")
        working.write(f">{fasta_transform_params.input_region.seq_id}\n")
        logger.debug(f"Fasta Fetch {fasta_transform_params.input_region}")

        sequence = fasta_file.fetch(fasta_transform_params.input_region.seq_id, fasta_transform_params.input_region.start - 1, fasta_transform_params.input_region.end)
        if len(sequence) < 50:
            logger.debug(f"Fasta Fetch = {sequence}")
        else:
            logger.debug(f"Fasta Fetch = {sequence[:25]} ... {sequence[-25:]}")
        g2g_utils.write_sequence(sequence, working)
        working.close()

        if fasta_transform_params.patch:
            logger.info("########################################################")
            logger.info("#################### PATCHING ##########################")
            logger.info("########################################################")

            p = fasta_patch.FastaPatchParams()
            p.input_region = fasta_transform_params.input_region

            # original input fasta file
            p.input_file = tmp_fasta

            # temporary directory
            p.temp_dir = fasta_transform_params.temp_dir

            p.output_file = tmp_fasta
            p.output_region = fasta_transform_params.input_region
            f_desc = (
                f"{fasta_transform_params.input_region.seq_id} "
                f"{fasta_transform_params.input_region.seq_id}:"
                f"{fasta_transform_params.input_region.start}-"
                f"{fasta_transform_params.input_region.end}"
            )
            p.output_header = fasta.FastaHeader(
                str(fasta_transform_params.input_region.seq_id),
                f_desc
            )

            # vci information
            p.vci_file = fasta_transform_params.vci_file
            p.vci_query = fasta_transform_params.vci_query
            p.reverse = fasta_transform_params.reverse

            p.gen_temp = False

            # offset
            p.offset = fasta_transform_params.offset

            p1 = fasta_patch.process_piece(p)
            logger.debug(f"\n\n\npatch output = {p1.output_file}\n\n\n")
            fasta_file = fasta.FastaFile(p1.output_file)
            fasta_transform_result.snp_count = p1.count
        else:
            fasta_file = fasta.FastaFile(tmp_fasta)

        offset = fasta_transform_params.offset
        reverse = fasta_transform_params.reverse

        logger.info(f"Transforming {fasta_transform_params.input_region}...")
        region = g2g.parse_region(fasta_transform_params.vci_query)
        logger.info(f"Finding VCI mappings for {region}")

        logger.info("PARSING VCI FILE...")
        vci_file = vci.VCIFile(fasta_transform_params.vci_file, seq_ids=[region.seq_id])
        vci_file.parse(reverse=reverse)
        logger.info("DONE PARSING!")

        mappings = vci_file.find_mappings(region.seq_id, fasta_transform_params.output_region.start-1, fasta_transform_params.output_region.end)
        fasta_out = open(fasta_transform_params.output_file, "w")

        if mappings is None:
            logger.info("This region was deleted")
            logger.info("TODO: dump the fasta sequence here")
            logger.info("fasta_transform_params.output_region.seq_id={}".format(fasta_transform_params.output_region.seq_id))
            logger.info("fasta_transform_params.output_file={}".format(fasta_transform_params.output_file))

            if fasta_transform_params.full_file:
                out_header = ">{} {}:{}-{} from|{}:{}-{}\n".format(fasta_transform_params.output_region.seq_id, fasta_transform_params.input_region.seq_id, fasta_transform_params.input_region.start, fasta_transform_params.input_region.end, fasta_transform_params.input_region.seq_id, fasta_transform_params.input_region.start, fasta_transform_params.input_region.end)
            else:
                out_header = ">{}:{}-{} from|{}:{}-{}\n".format(fasta_transform_params.output_region.seq_id, fasta_transform_params.input_region.start, fasta_transform_params.input_region.end, fasta_transform_params.input_region.seq_id, fasta_transform_params.input_region.start, fasta_transform_params.input_region.end)

            fasta_out.write(out_header)
            partial_seq = fasta_file.fetch(fasta_transform_params.input_region.seq_id, fasta_transform_params.input_region.start-1, fasta_transform_params.input_region.end)
            if len(partial_seq) < 50:
                logger.debug("Fasta Fetch = {}".format(partial_seq))
            else:
                logger.debug("Fasta Fetch = {} ... {}".format(partial_seq[:25], partial_seq[-25:]))
            g2g_utils.write_sequence(partial_seq, fasta_out)

            return fasta_transform_result

        first = True
        fasta_transform_result.ins_count = 0
        fasta_transform_result.del_count = 0

        new_start_pos = mappings[0].to_start
        new_end_pos = mappings[-1].to_end

        logger.debug("new_start_pos={}".format(new_start_pos))

        last_pos = 0
        new_sequence = StringIO()

        #logger.debug("index of '>' is {}".format(new_sequence.getvalue().find(">")))

        start = mappings[0].from_start
        logger.debug("Setting start to {} (mappings[0].from_start)".format(start))

        found = False
        new_sequence_len = 0

        logger.debug("start={}".format(start))
        logger.debug("last_pos={}".format(last_pos))
        logger.debug("VCI Fetch {}".format(fasta_transform_params.vci_query))

        for line in vci_file.fetch(fasta_transform_params.vci_query, parser=pysam.asTuple()):
            if line[5] == ".":
                continue

            found = True

            logger.debug("LINE: {}".format(line))
            #new_sequence_value = new_sequence.getvalue()

            #if len(new_sequence_value) > 50:
            #    logger.debug("current={}...{}".format(new_sequence_value[:25], new_sequence_value[-25:]))
            #else:
            #    logger.debug("current={}".format(new_sequence_value))

            # chromosome, position, shared_bases, deleted_bases, inserted_bases, fragment_size

            shared_bases = line[2]
            deleted_bases = line[3 if not reverse else 4]
            deleted_bases_length = len(deleted_bases) if deleted_bases != "." else 0
            inserted_bases = line[4 if not reverse else 3]
            inserted_bases_length = len(inserted_bases) if inserted_bases != "." else 0
            fragment_size = int(line[5])

            if first:
                # fragment_size = (current_pos + shared_bases_length) -
                #                 (previous_pos + previous_shared_length + previous_inserted_length)
                logger.debug("First result in query...")

                logger.debug("Adjusting last_pos from {} to {}".format(last_pos, start))
                last_pos = start

                logger.debug("Adjusting fragment_size from {} to {}".format(fragment_size, (int(line[1]) + len(shared_bases)) - (last_pos + 1 + 0)))
                fragment_size = (int(line[1]) + len(shared_bases)) - (last_pos + 1 + 0)

                first = False

                if fragment_size < 0:
                    continue

            logger.debug("last_pos={}".format(last_pos))
            logger.debug("offset={}".format(offset))
            logger.debug("fragment_size={}".format(fragment_size))
            #logger.debug("Fasta Fetch {}:{}-{} (0-based)".format(fasta_transform_params.input_region.seq_id, last_pos - offset, last_pos + fragment_size - offset))

            logger.debug("extracting... {}-{}".format(last_pos - offset, last_pos + fragment_size - offset))
            partial_seq = fasta_file.fetch(fasta_transform_params.input_region.seq_id, last_pos - offset, last_pos + fragment_size - offset)

            if len(partial_seq) < 50:
                logger.debug("Fasta Fetch = {}".format(partial_seq))
            else:
                logger.debug("Fasta Fetch = {} ... {}".format(partial_seq[:25], partial_seq[-25:]))

            new_sequence.write(partial_seq)
            new_sequence_len += fragment_size
            #logger.debug("index of '>' is {}".format(new_sequence.getvalue().find(">")))

            if inserted_bases_length > 0:
                # insertion
                logger.debug("INSERTION")
                new_sequence.write(inserted_bases)
                #logger.debug("index of '>' is {}".format(new_sequence.getvalue().find(">")))
                #logger.debug("{0}:{1}-{2} (Length: {3})".format(location.seqid, last_pos, last_pos + fragment_size, len(partial_seq)))
                #if len(partial_seq) > 100:
                #    logger.debug("{0}...{1}".format(partial_seq[:10], partial_seq[-10:]))
                #else:
                #    logger.debug(partial_seq)
                logger.debug("Adding {0}".format(inserted_bases))
                #logger.debug("SAME={0}, {1}".format(shared_bases, partial_seq[-(len(shared_bases)):]))

                fasta_transform_result.ins_count += inserted_bases_length
                new_sequence_len += inserted_bases_length

            if deleted_bases_length > 0:
                # deletion
                logger.debug("DELETION")
                last_pos += deleted_bases_length
                #logger.debug("skipping ahead {0} bases".format(deleted_bases_length))
                fasta_transform_result.del_count += deleted_bases_length

            #logger.debug("last_pos incremented by fragment_size, {} to {}".format(last_pos, last_pos + fragment_size))
            last_pos += fragment_size

            #logger.debug("LAST_POS={0}, INSERTIONS={1}, DELETIONS={2}, DIFF={3}".format(last_pos, fasta_transform_result.ins_count, fasta_transform_result.del_count, (fasta_transform_result.ins_count - fasta_transform_result.del_count)))

        if found:
            logger.debug("Fetching last bit of sequence")
            if last_pos >= fasta_transform_params.input_region.end:
                logger.debug("Nothing to fetch.. done")
            else:
                logger.debug("Fasta Fetch {}:{}-{} (0-based)".format(fasta_transform_params.input_region.seq_id, last_pos - offset, fasta_transform_params.input_region.end - offset))
                partial_seq = fasta_file.fetch(fasta_transform_params.input_region.seq_id, last_pos - offset, fasta_transform_params.input_region.end - offset)
                if len(partial_seq) < 50:
                    logger.debug("Fasta Fetch = {}".format(partial_seq))
                else:
                    logger.debug("Fasta Fetch = {} ... {}".format(partial_seq[:25], partial_seq[-25:]))
                new_sequence.write(partial_seq)
                #logger.debug("index of '>' is {}".format(new_sequence.getvalue().find(">")))
                new_sequence_len += (fasta_transform_params.input_region.end - last_pos)
        else:
            logger.debug("NO INDELS FOUND IN REGION")
            logger.debug("Fetching ONLY bit of sequence")
            logger.debug("start={}".format(start))
            logger.debug("offset={}".format(offset))
            logger.debug("fasta_transform_params.input_region.end={}".format(fasta_transform_params.input_region.end))

            logger.debug("Fasta Fetch {}:{}-{} (0-based)".format(fasta_transform_params.input_region.seq_id, start - offset, fasta_transform_params.input_region.end - offset))
            partial_seq = fasta_file.fetch(fasta_transform_params.input_region.seq_id, start - offset, fasta_transform_params.input_region.end - offset)
            #if len(partial_seq) < 50:
            #    logger.debug("Fasta Fetch = {}".format(partial_seq))
            #else:
            #    logger.debug("Fasta Fetch = {} ... {}".format(partial_seq[:25], partial_seq[-25:]))
            new_sequence.write(partial_seq)
            #logger.debug("index of '>' is {}".format(new_sequence.getvalue().find(">")))
            new_sequence_len += len(partial_seq) #(fasta_transform_params.input_region.end - offset) - 1

        logger.debug("the new_sequence_len={}".format(new_sequence_len))

        if fasta_transform_params.full_file:
            logger.debug("FULL FILE")
            out_header = ">{} {}:{}-{} from|{}:{}-{}\n".format(fasta_transform_params.output_region.seq_id, fasta_transform_params.input_region.seq_id, new_start_pos+1, new_start_pos + new_sequence_len, fasta_transform_params.input_region.seq_id, fasta_transform_params.input_region.start, fasta_transform_params.input_region.end)
        else:
            logger.debug("NOT FULL FILE")
            out_header = ">{}:{}-{} from|{}:{}-{}\n".format(fasta_transform_params.output_region.seq_id, new_start_pos+1, new_start_pos + new_sequence_len, fasta_transform_params.input_region.seq_id, fasta_transform_params.input_region.start, fasta_transform_params.input_region.end)

        logger.debug("WRITING HEADER: {}".format(out_header))
        fasta_out.write(out_header)
        #logger.debug("index of '>' is {}".format(new_sequence.getvalue().find(">")))
        g2g_utils.write_sequence(new_sequence.getvalue(), fasta_out)
        fasta_out.close()
    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except G2GRegionError as le:
        logger.debug("Unable to parse location, {0}".format(le))
        raise le
    except G2GValueError as e:
        logger.debug("Unable to parse alocation, {0}".format(e))
        raise e
    except G2GFastaError as e:
        logger.debug("Unable to parse blocation, {0}".format(e))
        raise e
    except TypeError as e:
        #g2g_utils._show_error()
        logger.debug("Unable to parse clocation, {0}".format(e))
        raise e
    except Exception as e:
        #g2g_utils._show_error()
        logger.debug("Unable to parse dlocation, {0}".format(e))
        raise e

    logger.info("Transforming complete for {}".format(fasta_transform_params.input_region))

    g2g_utils.delete_index_files(tmp_fasta)
    g2g_utils.delete_file(tmp_fasta)

    return fasta_transform_result
    #logger.info("Execution complete: {0}".format(g2g_utils.format_time(start_time, time.time())))


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    return process_piece(*args)


def prepare_fasta_transform(filename_fasta, filename_output):
    """
    Initialize fasta_transform variables

    """
    filename_output = g2g_utils.check_file(filename_output, "w")
    output_file_dir = os.path.abspath(os.path.dirname(filename_output))

    new_filename_output = filename_output

    # let's figure out what our output name will be
    if new_filename_output.lower().endswith(".gz"):
        # strip off .gz
        new_filename_output = new_filename_output[:-3]

    if not filename_output.lower().endswith(".fa"):
        raise G2GValueError("Expecting output filename extension to be either '.fa.gz' or '.fa'")

    g2g_utils.delete_index_files(new_filename_output)

    return new_filename_output


def process(filename_fasta, filename_vci, regions, filename_output=None, bgzip=False, reverse=False, num_processes=None, also_patch=False, debug_level=0):
    """
    Patch a Fasta file by replacing the bases where the SNPs are located in the VCF file.

    :param filename_fasta: name of the input Fasta file
    :type filename_fasta: string
    :param filename_g2g: name of the G2G file
    :type filename_g2g: string
    :param filename_output: name of the output Fasta file
    :type filename_output: string
    :param bgzip: compress file in BGZIP format
    :type bgzip: boolean
    :param reverse: reverse the G2G file
    :type reverse: boolean
    :param num_processes: the number of processes to spawn
    :type num_processes: int
    :return: Nothing
    """
    start = time.time()
    dump_fasta = False
    temp_directory = g2g_utils.create_temp_dir("transform_", dir=".")
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
            filename_output = g2g_utils.check_file(filename_output, "w")

            if not regions:
                filename_output = prepare_fasta_transform(filename_fasta, filename_output)
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
            filename_output = g2g_utils.gen_file_name(extension="fa", append_time=False, output_dir=temp_directory)
            dump_fasta = True
            logger.debug(f"Temporary fasta file: {filename_output}")

        fasta_file = fasta.FastaFile(filename_fasta)
        vci_file = vci.VCIFile(filename_vci)

        if fasta_file.is_diploid() and vci_file.is_haploid():
            raise G2GFastaError("Haploid VCI file and diploid Fasta file combination is not currently supported for transform")

        full_file = True

        if regions:
            full_file = False
            if len(regions) > 5:
                regions_str = ", ".join(l for l in map(str, regions[:5]))
                logger.warn(f"Regions: {regions_str} (showing 1st five)")
            else:
                regions_str = ", ".join(l for l in map(str, regions))
                logger.warn(f"Regions: {regions_str}")

        else:
            regions = []
            for chrom in fasta_file.references:
                regions.append(g2g.Region(chrom, 1, fasta_file.get_reference_length(chrom)))

        all_params = []
        for region in regions:
            logger.debug(f"region={region}")
            logger.debug(f"region.original_base={region.original_base}")

            params = FastaTransformParams()
            params.input_region = region
            params.input_file = filename_fasta
            params.temp_dir = temp_directory
            params.vci_file = filename_vci
            params.reverse = reverse
            params.offset = 0 if region.start <= 1 else region.start - 1
            params.patch = also_patch
            params.full_file = full_file
            params.debug_level = debug_level

            if fasta_file.is_haploid() and vci_file.is_diploid():
                # logger.info("*** Experimental ***")
                # logger.info("*** HAPLOID FASTA and DIPLOID VCI ***")
                logger.info("Fasta file is haploid and VCI file is diploid")
                params.output_file = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(region)+"_L", extension="fa", output_dir=temp_directory, append_time=False)

                if full_file:
                    params.output_region = g2g.Region(f"{region.seq_id}_L", region.start, region.end)
                    params.output_header = fasta.FastaHeader(f"{region.seq_id}_L", f"{region.seq_id}:{region.start}-{region.end}")
                else:
                    params.output_region = g2g.Region(f"{region.seq_id}_L", region.start, region.end)
                    params.output_header = fasta.FastaHeader(f"{region.seq_id}_L", f"{region.seq_id}_L:{region.start}-{region.end}")

                params.vci_query = f"{region.seq_id}_L:{region.start}-{region.end}"
                all_params.append(params)

                params_r = copy.deepcopy(params)
                params_r.output_file = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(region)+"_R", extension="fa", output_dir=temp_directory, append_time=False)

                if full_file:
                    params_r.output_region = g2g.Region(f"{region.seq_id}_R", region.start, region.end)
                    params_r.output_header = fasta.FastaHeader(f"{region.seq_id}_R", f"{region.seq_id}:{region.start}-{region.end}")
                else:
                    params_r.output_region = g2g.Region(f"{region.seq_id}_R", region.start, region.end)
                    params_r.output_header = fasta.FastaHeader(f"{region.seq_id}_R", f"{region.seq_id}_R:{region.start}-{region.end}")

                params_r.vci_query = f"{region.seq_id}_R:{region.start}-{region.end}"
                all_params.append(params_r)
            else:
                fmt_string = "HAPLOID" if vci_file.is_haploid() else "DIPLOID"
                logger.debug(f"VCI file and Fasta file are both: {fmt_string}")
                params.output_file = g2g_utils.gen_file_name(prefix=g2g_utils.location_to_filestring(region), extension="fa", output_dir=temp_directory, append_time=False)

                if full_file:
                    params.output_region = g2g.Region(region.seq_id, region.start, region.end)
                    params.output_header = fasta.FastaHeader(region.seq_id, f"{region.seq_id}:{region.start}-{region.end}")
                else:
                    params.output_region = g2g.Region(region.seq_id, region.start, region.end)
                    params.output_header = fasta.FastaHeader(f"{region.seq_id} {region.seq_id}:{region.start}-{region.end}", None)

                params.vci_query = f"{region.seq_id}:{region.start}-{region.end}"
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
            logger.warn(f"Processed {total_snp:,} SNPs total")

        logger.warn(f"Processed {total_ins:,} insertions total")
        logger.warn(f"Processed {total_del:,} deletions total")

        logger.debug("Temp files created, copy to master temp file and delete")
        g2g_utils.concatenate_files(files, filename_output, False)

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
                    g2g_utils.bgzip_and_index_file(filename_output, f"{filename_output}.gz", delete_original=True, file_format="fa")
                elif filename_output.lower().endswith(".gz"):
                    g2g_utils.bgzip_and_index_file(filename_output, filename_output, delete_original=True, file_format="fa")

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        raise G2GError(str(e))
    finally:
        # clean up the temporary files
        g2g_utils.delete_dir(temp_directory)
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warn(f"Transform complete: {fmt_time}")


#
# A way to combine vcf files
#

# standard library imports
from collections import OrderedDict
import multiprocessing
import os
import time

# 3rd party library imports
import pysam

# local library imports
from .exceptions import G2GError
from .exceptions import G2GValueError
from .exceptions import G2GVCFError
from .exceptions import KeyboardInterruptError
from . import fasta
from . import g2g
from . import g2g_utils
from . import vcf


class VCFFileInformation:
    def __init__(self, file_name=None, discard_file=None, sample_index=None):
        self.file_name = file_name
        self.discard_file = discard_file
        self.sample_index = sample_index
        self.lines = 0

    def __str__(self):
        return f"{self.file_name}, {self.sample_index}"


class VCF2VCIInfo(object):
    def __init__(self):
        self.chromosome = None

        # list of VCFFileInformation
        self.vcf_files = []
        self.fasta_file = None

        # indel specific

        self.prev_next_ref_pos_right = 1
        self.prev_next_ref_pos_left = 1

        self.diploid = False
        self.passed = False
        self.quality = False
        self.vcf_keep = False

        self.output_file_left = None
        self.output_file_right = None

        self.stats_left = {}
        self.stats_right = {}
        self.debug_level = 0


def walk_vcfs_together(readers, **kwargs):

    if "vcf_record_sort_key" in kwargs:
        get_key = kwargs["vcf_record_sort_key"]
    else:
        get_key = lambda r: (r.contig, r.pos)

    nexts = []

    for reader in readers:
        try:
            if reader:
               nexts.append(reader.next())
            else:
                nexts.append(None)
        except StopIteration:
            nexts.append(None)

    min_k = (None,)
    while any([r is not None for r in nexts]):
        next_idx_to_k = dict((i, get_key(r)) for i, r in enumerate(nexts) if r is not None)
        keys_with_prev_contig = [k for k in next_idx_to_k.values() if k[0] == min_k[0]]

        if any(keys_with_prev_contig):
            min_k = min(keys_with_prev_contig)
        else:
            min_k = min(next_idx_to_k.values())

        min_k_idxs = set([i for i, k in next_idx_to_k.items() if k == min_k])
        yield [nexts[i] if i in min_k_idxs else None for i in range(len(nexts))]

        for i in min_k_idxs:
            try:
                nexts[i] = readers[i].next()
            except StopIteration:
                nexts[i] = None


#: TODO: utilize stats
def update_stats(stats, reason):
    if not stats:
        stats = {}

    if reason in stats:
        stats[reason] += 1
    else:
        stats[reason] = 1

    return stats


def process_piece(merge_info):
    stats = {}
    logger = g2g.get_logger(merge_info.debug_level)

    try:
        output_file_left = None
        output_file_right = None

        if merge_info.output_file_left:
            output_file_left = open(merge_info.output_file_left, "w")

        if merge_info.output_file_right:
            output_file_right = open(merge_info.output_file_right, "w")

        mi = ["L"]
        if merge_info.diploid:
            mi = ["L", "R"]

        logger.warn(f"Processing Chromosome {merge_info.chromosome}...")

        iterators = []
        discard_functions = []

        for i, file_info in enumerate(merge_info.vcf_files):
            vcf_tabix = pysam.TabixFile(file_info.file_name)
            try:
                vcf_iterator = vcf_tabix.fetch(
                    merge_info.chromosome, parser=pysam.asVCF()
                )
                iterators.append(vcf_iterator)
            except ValueError as ve:
                iterators.append(None)

            if file_info.discard_file:
                vcf_discard = open(file_info.discard_file, "w")

                def discard_record(rec):
                    vcf_discard.write(str(rec))
                    vcf_discard.write("\n")

                discard_functions.append(discard_record)
            else:
                discard_functions.append(lambda rec: None)

        n = 0

        line_numbers = 0
        for vcf_records in walk_vcfs_together(iterators):
            for i, vcf_record in enumerate(vcf_records):
                # logger.debug(vcf_record)
                if vcf_record is None:
                    continue
                # logger.debug(vcf_record.alt)
                # logger.debug(type(vcf_record.alt))

                gt = vcf.parse_gt_tuple(
                    vcf_record, merge_info.vcf_files[i].sample_index
                )
                # logger.debug(gt)

                line_numbers = line_numbers + 1
                if gt.is_snp:
                    # snp
                    if merge_info.passed and "PASS" not in vcf_record.filter:
                        discard_functions[i](vcf_record)

                    # logger.debug("Processing SNP {}".format(vcf_record))
                    n += 1

                    if merge_info.quality and gt.fi == "0":
                        discard_functions[i](vcf_record)
                    elif gt.left is None or gt.right is None:
                        discard_functions[i](vcf_record)
                    else:
                        if merge_info.diploid:
                            # 0 i sthe same as REF and do not need
                            if gt.gt_left != 0:
                                output_file_left.write("{}_L\t{}\t{}\t{}\t{}\t{}\n".format(merge_info.chromosome, vcf_record.pos+1, ".", vcf_record.ref, gt.left, "."))
                            if gt.gt_right != 0:
                                output_file_right.write("{}_R\t{}\t{}\t{}\t{}\t{}\n".format(merge_info.chromosome, vcf_record.pos+1, ".", vcf_record.ref, gt.right, "."))
                        else:
                            if gt.gt_left == gt.gt_right and gt.gt_left != 0:
                                # ignore heterozygotes 0/1, 1/0, only process 0/0 and 1/1
                                # logger.debug("ACCEPTED")

                                # logger.debug("pos {} : ref {}, left {}, right {}".format(vcf_snp.pos, vcf_snp.ref, gt.left, gt.right))
                                output_file_left.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(merge_info.chromosome, vcf_record.pos+1, ".", vcf_record.ref, gt.left, "."))
                else:
                    # indel
                    logger.debug(f"Processing INDEL {vcf_record}")

                    if merge_info.passed and "PASS" not in vcf_record.filter:

                        logger.debug("TOSSED: FILTERED ON PASS")
                        logger.debug(vcf_record)
                        stats = update_stats(stats, "FILTERED ON PASS")

                        discard_functions[i](vcf_record)
                        continue

                    elif merge_info.quality and gt.fi == "0":

                        # FI : Whether a sample was a Pass(1) or fail (0) based on FILTER values

                        logger.debug("TOSSED: FILTERED ON QUALITY")
                        logger.debug(vcf_record)
                        stats = update_stats(stats, "FILTERED ON QUALITY")
                        discard_functions[i](vcf_record)
                        continue

                    elif gt.left is None and gt.right is None:

                        logger.debug("TOSSED: NO STRAIN DATA")
                        logger.debug(vcf_record)
                        stats = update_stats(stats, "NO STRAIN DATA")
                        logger.debug(i)
                        logger.debug(type(vcf_record))
                        discard_functions[i](vcf_record)
                        continue

                    elif not merge_info.diploid and gt.left != gt.right:
                        # haploid or hexaploid
                        # gt must be equal

                        logger.debug("TOSSED: HETEROZYGOUS")
                        logger.debug(vcf_record)
                        stats = update_stats(stats, "HETEROZYGOUS")
                        discard_functions[i](vcf_record)
                        continue

                    # START L AND R, ONLY R IF DIPLOID

                    for l_or_r in mi:
                        #logger.debug("******************")
                        #logger.debug(l_or_r)
                        lr_out = ""
                        if l_or_r == "L":
                            #logger.debug("->LEFT")
                            lr_out = "_L" if merge_info.diploid else ""
                            alt_seq = str(gt.left)
                            stats = merge_info.stats_left
                            output_file = output_file_left
                            prev_next_ref_pos = merge_info.prev_next_ref_pos_left
                        else:
                            #logger.debug("->RIGHT")
                            lr_out = "_R" if merge_info.diploid else ""
                            alt_seq = str(gt.right)
                            stats = merge_info.stats_right
                            output_file = output_file_right
                            prev_next_ref_pos = merge_info.prev_next_ref_pos_right

                        logger.debug(f"prev_next_ref_pos={prev_next_ref_pos}")

                        if gt.ref == alt_seq:

                            logger.debug("TOSSED, REF AND ALT ARE EQUAL")
                            logger.debug(vcf_record)
                            stats = update_stats(stats, "REF AND ALT ARE EQUAL")
                            discard_functions[i](vcf_record)
                            continue

                        orig_alt_seq = alt_seq

                        logger.debug(f"SAMPLE: {vcf_record[merge_info.vcf_files[i].sample_index]}")
                        logger.debug(f"REF='{gt.ref}', ALT_L='{gt.left}', ALT_R='{gt.right}', POS={vcf_record.pos}")

                        position = vcf_record.pos + 1

                        ref_seq = str(gt.ref)
                        len_ref = len(ref_seq)
                        len_alt = len(alt_seq)

                        base_changes = len_ref - len_alt
                        base_pos_diff = 0

                        if position < prev_next_ref_pos:
                            logger.debug(f"TOSSED: VCF ROLLBACK: {vcf_record}")
                            logger.debug(vcf_record)

                            stats = update_stats(stats, "VCF ROLLBACK")
                            discard_functions[i](vcf_record)
                            continue

                        # find the position where the first base change is
                        for n in range(min(len_ref, len_alt)):
                            if ref_seq[n] != alt_seq[n]:
                                base_pos_diff = n
                                break

                        # if it is 0, take the minimum length
                        if base_pos_diff == 0:
                            base_pos_diff = min(len_ref, len_alt)

                        # add the base position difference
                        position += base_pos_diff

                        # recalculate the strings
                        shared_bases = ref_seq[:base_pos_diff]
                        ref_seq = ref_seq[base_pos_diff:]
                        alt_seq = alt_seq[base_pos_diff:]

                        dt = len(ref_seq)
                        dq = len(alt_seq)

                        next_ref_pos = position + len(ref_seq)
                        fragment_size = position - prev_next_ref_pos

                        """

                        logger.debug("           gt.ref: {0}".format(gt.ref))
                        logger.debug("          ref_seq: {0}".format(ref_seq))
                        logger.debug("               dt: {0}".format(dt))
                        logger.debug("           gt.alt: {0}".format(orig_alt_seq))
                        logger.debug("          alt_seq: {0}".format(alt_seq))
                        logger.debug("               dq: {0}".format(dq))
                        logger.debug("         position: {0}".format(position))
                        logger.debug("prev_next_ref_pos: {0}".format(prev_next_ref_pos))
                        logger.debug("     next_ref_pos: {0}".format(next_ref_pos))
                        logger.debug("    fragment_size: {0}".format(fragment_size))
                        logger.debug("     base_changes: {0}".format(base_changes))
                        logger.debug("    base_pos_diff: {0}".format(base_pos_diff))
                        logger.debug("     shared_bases: {0}".format(shared_bases))
                        """
                        # fix any 0 length
                        if fragment_size < 0:
                            # logger.debug("TOSSED: FRAGMENT: {0}".format(vcf_record))

                            stats = update_stats(stats, "FRAGMENT SIZE < 0")
                            discard_functions[i](vcf_record)
                            continue

                        if fragment_size != 0:
                            ref_str = ref_seq if ref_seq else "."
                            alt_str = alt_seq if alt_seq else "."
                            out = f"{merge_info.chromosome}{lr_out}\t{vcf_record.pos+1}\t{shared_bases}\t{ref_str}\t{alt_str}\t{fragment_size}\n"
                            logger.debug(out)
                            output_file.write(out)
                        else:
                            #
                            # THIS SHOULD NOT HAPPEN
                            #
                            raise G2GVCFError("Conflicting VCF entries")

                        stats = update_stats(stats, "ACCEPTED")

                        if l_or_r == "L":
                            merge_info.stats_left = stats
                            merge_info.prev_next_ref_pos_left = next_ref_pos
                            logger.debug(f"setting merge_info.prev_next_ref_pos_left={merge_info.prev_next_ref_pos_left}")
                        else:
                            merge_info.stats_right = stats
                            merge_info.prev_next_ref_pos_right = next_ref_pos
                            logger.debug(f"setting merge_info.prev_next_ref_pos_right={merge_info.prev_next_ref_pos_right}")

        if merge_info.output_file_left:
            output_file_left.close()

        if merge_info.output_file_right:
            output_file_right.close()

    except KeyboardInterrupt:
        raise KeyboardInterruptError()
    except Exception as e:
        g2g_utils.show_error()
        raise Exception("Unknown exception")

    return {
        "chrom": merge_info.chromosome,
        "stats": stats,
        "merge_info": merge_info,
        "line_numbers": line_numbers
    }


def wrapper(args):
    """
    Simple wrapper, useful for debugging.

    :param args: the arguments to process_piece
    :return: the same as process_piece
    """
    return process_piece(*args)


def create_vci_header(
        temp_directory: str,
        fasta_file: str,
        vcf_input_files: list[VCFFileInformation],
        strain: str,
        vcf_keep: bool,
        passed: bool,
        quality: bool,
        diploid: bool,
        num_processes: int
) -> str:
    """
    Create a VCI file header that contains meta information about the VCI file.

    Args:
        temp_directory: Directory to build the file.
        fasta_file: Fasta file used in VCI file creation.
        vcf_input_files: List of VCIFileInformation.
        strain: The strain used to make the file.
        vcf_keep: True to place troubling VCF lines in extra file.
        passed: True uses only VCF lines that have a PASS for the filter.
        quality: True to filter on quality, FI=PASS.
        diploid: Create diploid VCI file.
        num_processes: Specify the number of processes.

    """
    file = g2g_utils.gen_file_name(
        "header", output_dir=temp_directory, extension="", append_time=False
    )

    with open(file, "w") as fd:
        create_time = time.strftime("%m/%d/%Y %H:%M:%S")
        fd.write(f"##CREATION_TIME={create_time}\n")

        for vcf_file in vcf_input_files:
            fd.write(f"##INPUT_VCF={vcf_file.file_name}\n")

        fd.write(f"##FASTA_FILE={fasta_file}\n")
        fd.write(f"##STRAIN={strain}\n")
        fd.write(f"##VCF_KEEP={vcf_keep}\n")
        fd.write(f"##FILTER_PASSED={passed}\n")
        fd.write(f"##FILTER_QUALITY={quality}\n")
        fd.write(f"##DIPLOID={diploid}\n")
        fd.write(f"##PROCESSES={num_processes}\n")

        fasta_file = fasta.FastaFile(fasta_file)

        for c in fasta_file.references:
            fd.write(f"##CONTIG={c}:{fasta_file.get_reference_length(c)}\n")

        fd.write("#CHROM\tPOS\tANCHOR\tINS\tDEL\tFRAG\n")
    return file


def process(
        vcf_files: list[str],
        fasta_file: str,
        output_file: str,
        strain: str,
        vcf_keep: bool | None = False,
        passed: bool | None = False,
        quality: bool | None = False,
        diploid: bool | None = False,
        num_processes: int = None,
        bgzip: bool | None = False,
        debug_level: int | None = 0
) -> None:
    """
    Parse the VCF file and create a VCI file.

    Args
        vcf_files: Name of the VCF files.
        fasta_file: Name of the Fasta file associated with the VCF file.
        bed_file_out: Name of output file.
        strain: Which strain to process.
        vcf_keep: True to place troubling VCF lines in extra file.
        passed: True uses only VCF lines that have a PASS for the filter.
        quality: True to filter on quality, FI=PASS.
        diploid: Create diploid VCI file.
        num_processes: Specify the number of processes.
        bgzip: True tp bgzip the VCI file.
        debug_level: Debug level (0=WARN,1=INFO,2+=DEBUG).
    """
    start = time.time()
    logger = g2g.get_logger(debug_level)

    output_file = g2g_utils.check_file(output_file, "w")
    output_file_dir = os.path.dirname(output_file)

    vcf_file_inputs = []

    if vcf_files:
        for file_name in vcf_files:
            vcf_file = g2g_utils.check_file(file_name)
            logger.warn(f"VCF file: {vcf_file}")
            logger.debug("Checking for index file, creating if needed...")
            g2g_utils.index_file(
                original_file=vcf_file, file_format="vcf", overwrite=False
            )

            vcf_discard_file = None
            if vcf_keep:
                vcf_discard_file = os.path.join(
                    output_file_dir,
                    f"{os.path.basename(output_file)}.errors.vcf"
                )
                logger.warn(f"VCF indel discard file: {vcf_discard_file}")

            vcf_file_inputs.append(
                VCFFileInformation(vcf_file, vcf_discard_file)
            )

    if len(vcf_file_inputs) == 0:
        raise G2GValueError("No VCF files.")

    if not fasta_file:
        raise G2GValueError("No fasta file was specified.")

    if not strain:
        raise G2GValueError("No strain was specified.")

    if not num_processes:
        num_processes = multiprocessing.cpu_count()
    else:
        if num_processes <= 0:
            num_processes = 1

    logger.warn(f"Fasta File: {fasta_file}")
    logger.warn(f"Strain: {strain}")
    logger.warn(f"Pass filter on: {passed}")
    logger.warn(f"Quality filter on: {quality}")
    logger.warn(f"Diploid: {diploid}")
    logger.warn(f"Number of processes: {num_processes}")
    logger.warn(f"Output VCI File: {output_file}")

    # not all chromosomes/seq_id will be processed if not in vcf file
    processed_seq_ids = {}

    temp_directory = g2g_utils.create_temp_dir("vcf2vci", dir=".")
    logger.debug(f"Temp directory: {temp_directory}")

    header_file = create_vci_header(
        temp_directory, fasta_file, vcf_file_inputs, strain,
        vcf_keep, passed, quality, diploid, num_processes
    )

    for i, vcf_file in enumerate(vcf_file_inputs):
        tb_file = pysam.TabixFile(vcf_file.file_name)
        for h in tb_file.header:
            h = g2g_utils.s(h)
            if h[:6] == "#CHROM":
                try:
                    elems = h.split("\t")
                    samples = elems[9:]
                    samples = dict(zip(samples, (x for x in range(len(samples)))))
                    vcf_file_inputs[i].sample_index = samples[strain]
                except KeyError:
                    elems = h.split("\t")
                    valid_strains = ", ".join(elems[9:])
                    raise G2GVCFError(
                        f"Unknown strain '{strain}', "
                        f"valid strains are: {valid_strains}"
                    )

        for seq_id in tb_file.contigs:
            processed_seq_ids[seq_id] = False

    tmp_processed_seq_ids = OrderedDict()

    for k in g2g_utils.natsorted(processed_seq_ids.keys()):
        tmp_processed_seq_ids[k] = False

    processed_seq_ids = tmp_processed_seq_ids

    all_merge_info = []
    pool = None

    try:
        for c in processed_seq_ids:
            merge_info = VCF2VCIInfo()
            merge_info.chromosome = c

            merge_info.vcf_files = vcf_file_inputs
            merge_info.fasta_file = fasta_file

            merge_info.diploid = diploid
            merge_info.passed = passed
            merge_info.quality = quality
            merge_info.vcf_keep = vcf_keep
            merge_info.debug_level = debug_level

            if diploid:
                merge_info.output_file_left = g2g_utils.gen_file_name(
                    f"chr{c}.left", output_dir=temp_directory,
                    extension="vci", append_time=False
                )

                merge_info.output_file_right = g2g_utils.gen_file_name(
                    f"chr{c}.right", output_dir=temp_directory,
                    extension="vci", append_time=False
                )

                g2g_utils.delete_file(merge_info.output_file_left)
                g2g_utils.delete_file(merge_info.output_file_right)
            else:
                merge_info.output_file_left = g2g_utils.gen_file_name(
                    f"chr{c}.right", output_dir=temp_directory,
                    extension="vci", append_time=False
                )

                g2g_utils.delete_file(merge_info.output_file_left)

            all_merge_info.append(merge_info)

        logger.info("Parsing VCF files...")

        args = zip(all_merge_info)
        pool = multiprocessing.Pool(num_processes)
        results = pool.map(wrapper, args)

        # parse results
        total = 0
        for r in results:
            total += r["line_numbers"]

        logger.debug("Combining temp files...")
        logger.info("Finalizing VCI file...")

        files = [header_file]
        if diploid:
            for mi in all_merge_info:
                files.append(mi.output_file_left)
                files.append(mi.output_file_right)

            g2g_utils.concatenate_files(files, output_file, True)

            if bgzip:
                g2g_utils.bgzip_and_index_file(
                    output_file, f"{output_file}.gz",
                    delete_original=True, file_format="vcf"
                )
        else:
            for mi in all_merge_info:
                files.append(mi.output_file_left)

            g2g_utils.concatenate_files(files, output_file, True)

            if bgzip:
                g2g_utils.bgzip_and_index_file(
                    output_file, output_file + ".gz",
                    delete_original=True, file_format="vcf"
                )

        # TODO: make sure stats are good and show statistics
        logger.warn("Parsed {0:,} total lines".format(total))
    except KeyboardInterruptError:
        pool.terminate()
        raise G2GError("Keyboard quit consumed")
    except KeyboardInterrupt:
        pool.terminate()
        raise G2GError("Execution halted")
    except Exception as e:
        g2g_utils.show_error()
        raise G2GError("Execution halted unknown error")
    finally:
        g2g_utils.delete_dir(temp_directory)
        fmt_time = g2g_utils.format_time(start, time.time())
        logger.warn(f"VCI creation complete: {fmt_time}")


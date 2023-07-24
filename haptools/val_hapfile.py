from __future__ import annotations

import os
from haptools import logging

from re import search
from pathlib import Path

from haptools import logging


LOGGER_NAME = "valhap"
LTS_SPEC = "0.2.0"
TRAIL = "\n>>>"


def tmpex(expectation: object, received: object) -> str:
    return f"Expected: {expectation}\nReceived: {received}"


class Line:
    def __init__(self, content: str, number: int):
        self.content: str = content
        self.number: int = number

        self.columns: list[str] = content.split()
        self.count: int = len(self.columns)

    def is_empty(self) -> bool:
        return self.count == 0

    def __getitem__(self, index: int) -> str:
        return self.columns[index]

    def __str__(self) -> str:
        return self.content


class HapFileIO:
    def __init__(self, filename: Path, logger=None):
        self.filename = filename
        self.log = logger or logging.getLogger(LOGGER_NAME)

    def lines(self, sorted: bool = True) -> list[Line]:
        buffer = open(self.filename)

        content = [
            Line(line.strip(), i + 1) for i, line in enumerate(buffer.readlines())
        ]
        content = list(filter(lambda line: not line.is_empty(), content))

        buffer.close()

        if not sorted:
            meta_limit = next(
                idx for idx, line in enumerate(content) if not line[0].startswith("#")
            )
            content = [
                line
                for idx, line in enumerate(content)
                if (not line[0].startswith("#")) or idx < meta_limit
            ]

        # lol
        content.sort(key=lambda line: ord(line[0][0]))

        return content

    def validate_existence(self) -> bool:
        if not self.exists():
            self.log.error(f"The file {self.filename} does not exist.")
            return False

        is_ok = True

        if not self.is_regular():
            self.log.error(f"Cannot read {self.filename}: Is not a regular file.")
            is_ok = False

        if not self.is_readable():
            self.log.error(f"Cannot read {self.filename}: Insufficient permissions.")
            is_ok = False

        return is_ok

    def exists(self) -> bool:
        return os.path.exists(self.filename)

    def is_regular(self):
        return os.path.isfile(self.filename)

    def is_readable(self) -> bool:
        return os.access(self.filename, os.R_OK)


class HapFileValidator:
    errc : int = 0
    warc : int = 0

    # H CHROM START END ID
    MANDATORY_HAPLOTYPE_COLUMN_COUNT: int = 5

    # R CHROM START END ID LN
    MANDATORY_REPEAT_COLUMN_COUNT: int = 5

    # V CHROM START END ID CHROM LN
    MANDATORY_VARIANT_COLUMN_COUNT: int = 6

    # # version <version>
    MANDATORY_VERSION_COLUMNS: int = 3

    # #X Name Type [Description]
    MANDATORY_DEFINITION_COLUMNS = 3

    KEY_HAPLOTYPE: int = 0
    KEY_REPEAT: int = 1
    KEY_VARIANT: int = 2
    KEY_VARIANT_SRC: int = 9

    NAME_HAPLOTYPE = "Haplotype"
    NAME_REPEAT = "Repeat"
    NAME_VARIANT = "Variant"

    KEY_KEY: str = "HT::Key"
    KEY_CHROMOSOME: str = "HT::Chromosome"
    KEY_START: str = "HT::Start"
    KEY_END: str = "HT::End"
    KEY_ID: str = "HT::ID"
    KEY_ALLELE: str = "HT::Allele"

    #DEFAULT_HEADER: dict[int, dict[str, type]] = {
    DEFAULT_HEADER: dict = {
        KEY_HAPLOTYPE: {},
        KEY_REPEAT: {},
        KEY_VARIANT: {},
    }

    #EMPTY_TYPES: dict[int, list[type]] = {
    EMPTY_TYPES: dict = {
        KEY_HAPLOTYPE: [],
        KEY_REPEAT: [],
        KEY_VARIANT: [],
    }

    #EMPTY_DATA: dict[int, list[Line]] = {
    EMPTY_DATA: dict = {
        KEY_HAPLOTYPE: [],
        KEY_REPEAT: [],
        KEY_VARIANT: [],
    }

    #EMPTY_HRIDS: dict[int, dict[str, Line]] = {
    EMPTY_HRIDS: dict = {
        KEY_HAPLOTYPE: {},
        KEY_REPEAT: {},
    }

    #EMPTY_VRIDS: dict[str, dict[str, Line]] = {}
    EMPTY_VRIDS: dict = {}

    #EMPTY_META: list[Line] = []
    EMPTY_META: list = []

    def __init__(self, logger=None):
        self.log = logger or logging.getLogger(LOGGER_NAME)

        self.vars_ex: dict[int, dict[str, type]] = HapFileValidator.DEFAULT_HEADER
        self.types_ex: dict[int, list[type]] = HapFileValidator.EMPTY_TYPES

        self.meta: list[Line] = HapFileValidator.EMPTY_META
        self.data: dict[int, list[Line]] = HapFileValidator.EMPTY_DATA

        self.hrids: dict[int, dict[str, Line]] = HapFileValidator.EMPTY_HRIDS
        self.vrids: dict[str, dict[str, Line]] = HapFileValidator.EMPTY_VRIDS

        self.referenced_chromosomes: set[str] = set()

    def extract_and_store_content(self, file: HapFileIO, sorted: bool = True):
        lines = file.lines(sorted=sorted)

        self.extract_meta_lines(lines)
        self.extract_data_lines(lines)

    def extract_meta_lines(self, lines: list[Line]):
        header_limit = next(
            i for i, line in enumerate(lines) if not line[0].startswith("#")
        )
        self.meta_lines = lines[:header_limit]

    def extract_data_lines(self, lines: list[Line]):
        limits = [0, 0, 0]
        for i, char in enumerate(["H", "R", "V"]):
            limits[i] = next(
                i for i, line in enumerate(lines) if line[0].startswith(char)
            )

        ln = [
            lines[limits[0] : limits[1]],
            lines[limits[1] : limits[2]],
            lines[limits[2] :],
        ]

        for i in range(HapFileValidator.KEY_HAPLOTYPE, HapFileValidator.KEY_VARIANT + 1):
            self.data[i] = ln[i]

    #
    # Version Validation
    #

    def validate_version_declarations(self):
        versions = self.extract_version_declarations()
        if len(versions) == 0:
            self.log.warn(
                f"{TRAIL} No version declaration found. Assuming to use the latest"
                " version."
            )
            HapFileValidator.warc += 1
            return

        self.validate_version_format(versions[-1])

    def extract_version_declarations(self) -> list[Line]:
        decls = list(
            filter(lambda x: x.count > 1 and x[1] == "version", self.meta_lines)
        )

        if len(decls) > 1:
            self.log.warn(
                f"{TRAIL} Found more than one version declaration. Using the last"
                " instance. Each is its own warning."
            )

            for decl in decls:
                HapFileValidator.warc += 1
                self.lwfl("", decl, sep="")

        return decls

    def validate_version_format(self, version: Line):
        if version.count < 3:
            self.leexfl(
                "Not enough columns in version declaration",
                HapFileValidator.MANDATORY_DEFINITION_COLUMNS,
                version.count,
                version,
            )
            self.warnskip(version)

            HapFileValidator.errc += 1
            return

        if search(r"\d+\.\d+\.\d+", version[2]) == None:
            self.leexfl(
                "Version is incorrectly formatted",
                "'x.x.x' where 'x' is an integer",
                version[2],
                version,
            )

            HapFileValidator.errc += 1

    #
    # Column additions
    #

    def validate_column_additions(self):
        additions = self.find_column_additions()

        for i, k in enumerate(["#H", "#R", "#V"]):
            self.add_column_additions_to_header(
                i, list(filter(lambda line: line[0] == k, additions))
            )

    def find_column_additions(self) -> list[Line]:
        additions = list(
            filter(lambda line: search(r"#[H|R|V]", line[0]) != None, self.meta_lines)
        )

        invalid_lines = [
            x for x in self.meta_lines if x not in additions and len(x[0]) > 1
        ]

        for ln in invalid_lines:
            self.leexfl(
                "Invalid column addition type.",
                "A column addition for 'H', 'R', or 'V'",
                f"A column addition for '{ln[0][1]}', whose type doesn't exist",
                ln,
            )
            HapFileValidator.errc += 1

        return additions

    def add_column_additions_to_header(self, tp: int, additions: list[Line]):
        for addition in additions:
            if addition.count < 3:
                self.leexfl(
                    "Insufficient columns for extra column definition",
                    HapFileValidator.MANDATORY_DEFINITION_COLUMNS,
                    addition.count,
                    addition,
                )
                self.warnskip(addition)

                HapFileValidator.errc += 1
                return

            ptp = self.retrieve_column_addition_data_type(addition)

            if ptp == object:
                self.warnskip(addition)
                continue

            self.vars_ex[tp].update({addition[1]: ptp})
            self.types_ex[tp].append(ptp)

    def retrieve_column_addition_data_type(self, addition: Line) -> type:
        tp = addition[2]

        if tp == "d":
            return int

        if search(r"\.\d+f", addition.content) != None:
            return float

        if tp == "s":
            return str

        self.leexfl(
            "Could not parse type for column addition",
            "One of: 'd', 's', '.xf' (where 'x' is an integer)",
            f"{addition[2]}",
            addition,
        )
        
        HapFileValidator.errc += 1
        return object

    #
    # Minimum Requirements
    #

    def validate_columns_fulfill_minreqs(self):
        self.validate_haplotypes()
        self.validate_repeats()
        self.validate_variants()

    def validate_haplotypes(self):
        for line in self.data[HapFileValidator.KEY_HAPLOTYPE]:
            has_min_cols = self.check_has_min_cols(
                line, HapFileValidator.MANDATORY_HAPLOTYPE_COLUMN_COUNT
            )

            self.check_start_and_end_positions(line)

            if not has_min_cols:
                self.lwexfl(
                    "Cannot check for variant references: Insufficient columns",
                    "A mandatory 5 columns for haplotyes",
                    line.count,
                    line,
                )
                self.warnskip(line)

                HapFileValidator.warc += 1
                return

            variant_refs = self.vrids.get(line[4])

            if variant_refs == None:
                self.leexfl(
                    f"Haplotype ID '{line[4]}' is not associated to any variants",
                    f"A variant association for Haplotype ID '{line[4]}'",
                    "No association",
                    line,
                )

                HapFileValidator.errc += 1
                return

    def validate_repeats(self):
        for line in self.data[HapFileValidator.KEY_REPEAT]:
            self.check_has_min_cols(line, HapFileValidator.MANDATORY_REPEAT_COLUMN_COUNT)

            self.check_start_and_end_positions(line)

    def validate_variants(self):
        for line in self.data[HapFileValidator.KEY_VARIANT]:
            self.check_has_min_cols(line, HapFileValidator.MANDATORY_VARIANT_COLUMN_COUNT)

            self.check_start_and_end_positions(line)
            self.check_variant_alleles(line)

    def check_has_min_cols(self, line: Line, min: int) -> bool:
        if line.count < min:
            self.leexfl(
                "Invalid amount of mandatory columns in definition.",
                f"At least {min}",
                line.count,
                line,
            )

            HapFileValidator.errc += 1
            return False

        return True

    def check_start_and_end_positions(self, line: Line):
        if line.count < 3:
            self.lwfl(
                "Cannot validate start and end positions: Insufficient columns", line
            )
            self.warnskip(line)

            HapFileValidator.warc += 1
            return

        f = False

        if not line[2].isdigit():
            self.leexfl(
                "Cannot convert start position to integer",
                "Integer values for the start position",
                line[2],
                line,
            )

            HapFileValidator.errc += 1
            f = True

        if not line[3].isdigit():
            self.leexfl(
                "Cannot convert end position to integer",
                "Integer values for the end position",
                line[3],
                line,
            )

            HapFileValidator.errc += 1
            f = True

        if f:
            self.lwfl(
                "Cannot test for correct position order due to previous errors"
                " (Inconvertible integers)",
                line,
            )
            self.warnskip(line)

            HapFileValidator.warc += 1
            return

        start = int(line[2])
        end = int(line[3])

        if start > end:
            self.leexfl(
                "Start position is greater than the end position",
                f"Start to be positioned at or before the end",
                f"{start} > {end} | Difference of {start - end}",
                line,
            )

            HapFileValidator.errc += 1

        if line.count < 5:
            self.lwexfl(
                "Cannot perform position validations against variant definitions:"
                " Insufficient columns.",
                5,
                line.count,
                line,
            )
            self.warnskip(line)

            HapFileValidator.warc += 1
            return

        variant_refs = self.vrids.get(line[4])

        if variant_refs == None:
            return

        for id, ln in variant_refs.items():
            if not ln[2].isdigit():
                self.leexfl(
                    "Variant start position cannot be converted to an integer.",
                    "An integer",
                    ln[2],
                    ln,
                )
                self.warnskip(line)

                HapFileValidator.errc += 1
                return

            if not ln[3].isdigit():
                self.leexfl(
                    "Variant end position cannot be converted to an integer.",
                    "An integer",
                    ln[3],
                    ln,
                )
                self.warnskip(line)

                HapFileValidator.errc += 1
                return

            vstart = int(ln[2])
            vend = int(ln[3])

            if vstart < start:
                self.leexfl(
                    "Variant start position cannot be prior to the start position of"
                    " its haplotype.",
                    "The variant to start after or when the haplotype does",
                    f"[Variant] {vstart} < [Haplotype] {start} | Difference of"
                    f" {start - vstart}",
                    line,
                )
                self.log.warn(f"At Line #{ln.number}: {ln}")

                HapFileValidator.errc += 1

            if vend > end:
                self.leexfl(
                    "Variant end position cannot be after than the end position of its"
                    " haplotype.",
                    "The variant to end before or when the haplotype does",
                    f"[Variant] {vend} > [Haplotype] {end} | Difference of"
                    f" {vend - end}",
                    line,
                )
                self.log.warn(f"At Line #{ln.number}: {ln}")

                HapFileValidator.errc += 1

    def check_variant_alleles(self, line: Line):
        if line.count < HapFileValidator.MANDATORY_VARIANT_COLUMN_COUNT:
            self.lwexfl(
                "Cannot test for variant allele type: Not enough columns.",
                HapFileValidator.MANDATORY_VARIANT_COLUMN_COUNT,
                line.count,
                line,
            )
            self.warnskip(line)

            HapFileValidator.warc += 1
            return

        if line[5].upper() not in ["A", "C", "G", "T"]:
            self.leexfl(
                "Invalid allele type in variant.",
                "One of 'A', 'C', 'G', 'T'",
                f"'{line[5]}'",
                line,
            )

            HapFileValidator.errc += 1

    #
    # ID Storage
    #

    def store_ids(self):
        for tp in range(2):
            for line in self.data[tp]:
                self.store_hrid(tp, line)

        for line in self.data[HapFileValidator.KEY_VARIANT]:
            self.store_variant_id(line)

    def store_hrid(self, tp: int, line: Line):
        should_skip = False
        if line.count < 2:
            self.lwexfl(
                "Cannot extract chromosome ID: Insufficient columns.",
                "At least 1 column",
                line.count,
                line,
            )

            HapFileValidator.warc += 1
            should_skip = True

        if line.count < 5:
            self.lwexfl(
                "Cannot extract ID: Insufficient columns.",
                f"At least 5 for ID extraction",
                line.count,
                line,
            )

            HapFileValidator.warc += 1
            should_skip = True

        if should_skip:
            self.warnskip(line)
            return

        self.referenced_chromosomes.add(line[1])

        if line[4] in self.hrids[tp]:
            self.leexfl("Duplicate ID.", "A unique ID", f"'{line[4]}'", line)
            self.log.warn(
                f"Originally defined at: line #{self.hrids[tp][line[4]].number}"
            )

            self.warnskip(line)

            HapFileValidator.errc += 1
            return

        if line[4] in self.referenced_chromosomes:
            self.lwfl(f"ID '{line[4]}' is already registered as a chromosome.", line)
            self.warnskip(line)

            HapFileValidator.warc += 1
            return

        self.hrids[tp].update({line[4]: line})

    def store_variant_id(self, line: Line):
        if line.count < 5:
            self.lwexfl(
                "Cannot extract ID: Insufficient columns.",
                f"At least 5 for ID extraction",
                line.count,
                line,
            )
            HapFileValidator.warc += 1

        if not line[1] in self.vrids.keys():
            self.vrids.update({line[1]: {}})

        if line[4] in self.vrids[line[1]].keys():
            self.leexfl(
                "Duplicate variant in for a same haplotype ID.",
                "A unique ID per haplotype",
                f"'{line[4]}'",
                line,
            )
            self.log.warn(
                f"Originally defined at: line #{self.vrids[line[1]][line[4]].number}"
            )
            self.warnskip(line)

            HapFileValidator.errc += 1
            return

        if line[4] in self.referenced_chromosomes:
            self.lwfl(f"ID '{line[4]}' is already registered as a chromosome.", line)
            self.warnskip(line)

            HapFileValidator.warc += 1
            return

        self.vrids[line[1]].update({line[4]: line})

    #
    # Variant Validation
    #

    def validate_variants_against_haplotypes(self):
        self.validate_variant_ids()

    def validate_variant_ids(self):
        for haplotype, ids in self.vrids.items():
            no_haplotype = False
            for id, line in ids.items():
                if haplotype not in self.hrids[HapFileValidator.KEY_HAPLOTYPE].keys():
                    self.lefl(
                        f"Cannot link variant '{id}' to non-exisent haplotype"
                        f" '{haplotype}'",
                        line,
                    )
                    no_haplotype = True

                    HapFileValidator.errc += 1
                    continue

            if no_haplotype:
                self.log.warn(
                    f"{TRAIL} Define haplotype '{haplotype}' or fix the variant"
                    " haplotype reference"
                )

    #
    # Extra field validation
    #

    def validate_extra_fields(self):
        for tp in range(HapFileValidator.KEY_HAPLOTYPE, HapFileValidator.KEY_VARIANT + 1):
            excol_count = len(self.types_ex[tp])
            lines = self.data[tp]

            for line in lines:
                rs = 5 if tp != HapFileValidator.KEY_VARIANT else 6
                extras = line.count - rs
                if extras != excol_count:
                    self.lwexfl(
                        "Invalid amount of extra columns in line.",
                        excol_count,
                        extras,
                        line,
                    )

                    if extras < 0:
                        self.lefl("There aren't even enough mandatory columns", line)
                        HapFileValidator.warc += 1

                    self.warnskip(line)

                    HapFileValidator.warc += 1
                    continue

                for ptp, col in zip(self.types_ex[tp], line.columns[rs:]):
                    conv = self.determine_if_is_convertible(col, ptp)

                    if not conv:
                        self.leexfl(
                            "Value in extra column is not convertible to the associated"
                            " type",
                            f"A value that can be converted to a(n) {str(ptp)[8:-2]}",
                            col,
                            line,
                        )

                        HapFileValidator.errc += 1

    def determine_if_is_convertible(self, what: str, tp: type) -> bool:
        if tp == int:
            return what.isdigit()

        if tp == float:
            return search(r"\d*\.?\d+$", what) != None

        return tp == str

    #
    # Extra field reordering
    #

    def reorder_extra_fields(self):
        reordering_metalns = list(
            filter(
                lambda line: line.count > 1 and search("order[H|R|V]", line[1]) != None,
                self.meta_lines,
            )
        )

        for i, c in enumerate(["H", "R", "V"]):
            relevant = list(filter(lambda line: line[1][5] == c, reordering_metalns))

            if len(relevant) == 0:
                continue

            if len(relevant) > 1:
                self.log.warn(
                    f"Found multiple order{c} definition lines. Using the last"
                    " available one."
                )
                HapFileValidator.warc += 1

            ln = relevant[-1]

            self.reorder_field_types(i, ln)

    def reorder_field_types(self, tp: int, line: Line):
        extpc = len(self.vars_ex[tp].keys())
        exclc = line.count - 2

        if extpc != exclc:
            self.leexfl(
                "Not enough columns in extra column reordering", extpc, exclc, line
            )
            self.warnskip(line)

            HapFileValidator.errc += 1
            return

        s = False
        for col in line.columns[2:]:
            if not col in self.vars_ex[tp]:
                self.lefl(f"{col} has not been defined as an extra colunm", line)
                HapFileValidator.errc += 1
                s = True

        if s:
            self.warnskip(line)
            return

        self.types_ex[tp].clear()
        for col in line.columns[2:]:
            self.types_ex[tp].append(self.vars_ex[tp][col])

    #
    # Logging
    #

    def lefl(self, msg: str, line: Line, sep: str = "\n"):
        self.log.error(f"{TRAIL} {msg}{sep}At line #{line.number}: {line}")

    def lwfl(self, msg: str, line: Line, sep: str = "\n"):
        self.log.warn(f"{TRAIL} {msg}{sep}At line #{line.number}: {line}")

    def lwexfl(self, msg: str, exp: object, rec: object, line: Line, sep: str = "\n"):
        self.log.warning(
            f"{TRAIL} {msg}{sep}{tmpex(exp, rec)}{sep}At line #{line.number}: {line}"
        )

    def leexfl(self, msg: str, exp: object, rec: object, line: Line, sep: str = "\n"):
        self.log.error(
            f"{TRAIL} {msg}{sep}{tmpex(exp, rec)}{sep}At line #{line.number}: {line}"
        )

    def warnskip(self, line: Line):
        self.log.warning(f"Skipping line #{line.number}")


def is_hapfile_valid(filename: Path, sorted=True, logger=None) -> bool:
    log = logger

    if log == None:
        log = logging.getLogger(LOGGER_NAME)


    file = HapFileIO(filename, logger=log)
    errc = 0

    is_readable = file.validate_existence()

    if not is_readable:
        return False

    hapfile = HapFileValidator(logger=log)
    hapfile.extract_and_store_content(file, sorted=sorted)

    hapfile.store_ids()

    hapfile.validate_column_additions()

    hapfile.validate_columns_fulfill_minreqs()
    hapfile.validate_variant_ids()

    hapfile.reorder_extra_fields()

    hapfile.validate_extra_fields()

    hapfile.validate_version_declarations()

    log.info(f"Completed HapFile validation with {HapFileValidator.errc} errors and {HapFileValidator.warc} warnings.")

    return HapFileValidator.errc == 0 and HapFileValidator.warc == 0



from pathlib import Path
import re
from sys import stderr, argv
from os import R_OK, access


def error_expected(what : str, got : str, linen : int) -> None:
    print(f">>> Expected {what} but got {got}", file = stderr)
    print(f">>> At line {linen}\n", file = stderr)


def error_expecting_cols(cols : list[str]) -> None:
    print(f">>>>> Expecting: {cols}\n", file = stderr)


def is_line_meta(line : str) -> bool:
    return line.startswith("#")


def is_line_regular(line : str) -> bool:
    return not is_line_meta(line)


def is_convertible_to_type(tp : type, what : str) -> bool:
    if tp == int:
        return what.isdigit()

    elif tp == float:
        return what.isdigit()

    elif tp == str:
        return True

    return False



def assert_file_exists(path : Path) -> None:
    if (path.exists()):
        return

    print(f">>> Could not open {path}: The file does not exist", file = stderr)


def assert_file_is_regular(path : Path) -> None:
    if path.is_file():
        return

    print(f">>> Failed to read {path}: It is not a regular file.")
    quit(1)


def assert_file_is_readable(path : Path) -> None:
    if access(path, R_OK):
        return

    print(f">>> Failed to read {path}:Not enough permissions.")
    quit(1)


def read_file_lines_protected(filename : str) -> list[str]:
    content : list[str]
    path : Path = Path(filename)

    assert_file_exists(path)
    assert_file_is_regular(path)
    assert_file_is_readable(path)

    buffer = open(path)

    content = buffer.readlines()

    buffer.close()

    return content


class Columns:


    def __init__(self, content : list[str]):
        self.content : list[str] = content
        self.count   : int       = len(content)


    def assert_length_eq(self, length : int) -> bool:
        if self.count == length:
            return True

        return False


    def assert_length_gte(self, length : int) -> bool:
        if self.count >= length:
            return True

        return False


    def get(self, index : int) -> str:
        if index >= self.count:
            return ""

        return self.content[index]


class Line:


    def __init__(self, number : int, content : str):
        self.number  : int            = number
        self.content : str            = content
        self.columns : Columns | None = None

        self.fatal   : int = 0
        self.warning : int = 0


    def split_and_save(self) -> None:
        self.columns = Columns(self.content.split())


    def as_columns(self) -> Columns:
        if (self.columns == None):
            self.split_and_save()
            assert self.columns != None

        return self.columns


    def is_flawed(self) -> bool:
        return self.is_wrong() or self.is_ill_formed()


    def is_wrong(self) -> bool:
        return self.fatal > 0


    def is_ill_formed(self) -> bool:
        return self.warning > 0


    def err(self) -> None:
        self.fatal += 1

    
    def warn(self) -> None:
        self.warning += 1

    def display_flaws_if_present(self) -> None:
        if not self.is_flawed():
            return

        print(f">>>>>>> {self.fatal} error(s) and {self.warning} warning(s) emmited for line #{self.number}", file = stderr)
        print(f">>>>>>> {self.content}", file = stderr)


class HapFile:


    KEYWORD_VERSION       = r"version"
    KEYWORD_VERSION_REGEX = r"\d+\.\d+\.\d+"
    KEYWORD_ORDER_REGEX   = r"order."

    NONE      : int = 0
    HAPLOTYPE : int = 1
    REPEAT    : int = 2
    VARIANT   : int = 3

    MANDATORY_COLUMNS_HAPLOTYPE   : int = 5
    MANDATORY_COLUMNS_REPEAT      : int = 5
    MANDATORY_COLUMNS_VARIANT     : int = 6
    MANDATORY_COLUMNS_FIELD_DEF   : int = 3
    MANDATORY_COLUMNS_VERSION_DEF : int = 3

    KEY_LINE_TYPE      = "HT_LINE_TYPE"
    KEY_CHROMOSOME_ID  = "HT_CHROMOSOME_ID"
    KEY_START_POSITION = "HT_START_POSITION"
    KEY_END_POSITION   = "HT_END_POSITION"
    KEY_ID             = "HT_ID"
    KEY_ALLELE         = "HT_ALLELE"


    CHARACTER_TYPE_ASSOCIATIONS : dict[str, int] = {
        "H" : HAPLOTYPE,
        "R" : REPEAT,
        "V" : VARIANT
    }

    CHARACTER_PYTYPE_ASSOCIATIONS : dict[str, type] = {
        "d" : int,
        "f" : float,
        "s" : str,
    }


    DEFAULT_HEADER : dict[int, list[tuple[str, type]]] = {
        HAPLOTYPE: [
            (KEY_LINE_TYPE,      str),
            (KEY_CHROMOSOME_ID,  str),
            (KEY_START_POSITION, int),
            (KEY_END_POSITION,   int),
            (KEY_ID,             str)
        ],

        REPEAT: [
            (KEY_LINE_TYPE,      str),
            (KEY_CHROMOSOME_ID,  str),
            (KEY_START_POSITION, int),
            (KEY_END_POSITION,   int),
            (KEY_ID,             str)
        ],

        VARIANT: [
            (KEY_LINE_TYPE,      str),
            (KEY_CHROMOSOME_ID,  str),
            (KEY_START_POSITION, int),
            (KEY_END_POSITION,   int),
            (KEY_ID,             str),
            (KEY_ALLELE,         str)
        ]
    }

    DEFAULT_TABLE : dict[int, list[list[str]]] = {
        HAPLOTYPE : [],
        REPEAT    : [],
        VARIANT   : []
    }

    @staticmethod
    def get_associated_hapfile_type_from_str(s : str) -> int | None:
        return HapFile.CHARACTER_TYPE_ASSOCIATIONS.get(s.upper())

    @staticmethod
    def get_associated_pytype_from_str(s : str) -> type | None:
        return HapFile.CHARACTER_PYTYPE_ASSOCIATIONS.get(s[len(s) - 1])


    def __init__(self):
        self.header  : dict[int, list[tuple[str, type]]] = HapFile.DEFAULT_HEADER
        self.tensor  : dict[int, list[list[str]]]        = HapFile.DEFAULT_TABLE
        self.version : str | None = None

        self.fatal_errors : int = 0
        self.warnings     : int = 0


    #
    # Reading
    #


    def read_file(self, filename : str) -> None:
        file_content = read_file_lines_protected(filename)

        header_lines : list[Line] = []
        data_lines   : list[Line] = []

        linen = 0

        for line in file_content:
            linen += 1

            if line.isspace(): continue

            if is_line_meta(line):
                header_lines.append(Line(linen, line))
            else:
                data_lines.append(Line(linen, line))


        self.read_into_header(header_lines)
        self.read_into_matrix(data_lines)


    #
    # Header
    #


    def read_into_header(self, values : list[Line]) -> None:
        for line in values:
            self.parse_meta_line(line)
            line.display_flaws_if_present()


    def parse_meta_line(self, line : Line):
        columns : Columns = line.as_columns()

        if len(columns.get(0)) < 2:
            self.parse_comment_or_meta(line)

        else:
            self.parse_column_addition(line)

    #
    # C1 Is Just "#"
    #


    def parse_comment_or_meta(self, line : Line) -> None:
        columns : Columns = line.as_columns()

        if columns.count < 2: return
        metatype_column = columns.get(1)

        if metatype_column == HapFile.KEYWORD_VERSION:
            self.parse_version(line)

        elif re.search(HapFile.KEYWORD_ORDER_REGEX, metatype_column):
            self.set_column_order(line)


    #
    # Version Key Present
    #


    def parse_version(self, line: Line) -> None:
        columns = line.as_columns()

        s = columns.assert_length_eq(HapFile.MANDATORY_COLUMNS_VERSION_DEF)
        if not s:
            line.err()
            error_expected(f"{HapFile.MANDATORY_COLUMNS_VERSION_DEF} columns for version definition", str(columns.count), line.number)
            error_expecting_cols(["Hash (#)", "version", "<x.x.x>"])
            self.skip_due_to_errors_emmited(line)
            return

        version = line.as_columns().get(2)
        if (re.search(HapFile.KEYWORD_VERSION_REGEX, version)) == None:
            error_expected("a version whose format conforms to \"x.x.x\" where \"x\" is an integer", f"\"{version}\"", line.number)
            line.warn()

        self.version = version

    #
    # orderX present
    #


    def set_column_order(self, line : Line) -> None:
        columns : Columns = line.as_columns()

        order_x          = columns.get(1)
        hapfile_type_str = order_x[5:]

        tp = HapFile.CHARACTER_TYPE_ASSOCIATIONS.get(hapfile_type_str)

        if (tp == None):
            error_expected(f"one of {list(HapFile.CHARACTER_TYPE_ASSOCIATIONS.keys())} as the type for an order definition", f"\"{hapfile_type_str}\"", line.number)

        print("CALL set_column_order <UNIMPLEMENTED>")


    #
    # C1 is #[...]
    #


    def parse_column_addition(self, line: Line) -> None:
        columns : Columns = line.as_columns()

        success = columns.assert_length_gte(HapFile.MANDATORY_COLUMNS_FIELD_DEF)
        if not success:
            line.err()
            error_expected(f"{HapFile.MANDATORY_COLUMNS_FIELD_DEF} columns for extra field definition", str(columns.count), line.number)
            error_expecting_cols(["Hash & Type (#X)", "Name", "Data Type Format (s, d, .xf)", "Optional Description"])
            self.skip_due_to_errors_emmited(line)
            return

        hapfile_type_str = columns.get(0)[1:]
        hapfile_type     = HapFile.get_associated_hapfile_type_from_str(hapfile_type_str)

        if (hapfile_type == None):
            line.warn()
            error_expected(f"one of {list(HapFile.CHARACTER_TYPE_ASSOCIATIONS.keys())} for the line type when adding an extra field", f"\"{hapfile_type_str}\"", line.number)

        python_type = HapFile.get_associated_pytype_from_str(columns.get(2))

        if python_type == None:
            line.warn()
            error_expected(f"one of {list(HapFile.CHARACTER_PYTYPE_ASSOCIATIONS.keys())} for the data type when adding an extra field", f"\"{columns.get(2)}\"", line.number)

        if (line.is_flawed()):
            self.skip_due_to_errors_emmited(line)
            return

        assert hapfile_type != None
        assert python_type  != None

        self.header[hapfile_type].append((columns.get(1), python_type))


    #
    # Matrix
    #

    
    def read_into_matrix(self, values : list[Line]) -> None:
        for line in values:
            self.parse_data_line(line)
            line.display_flaws_if_present()


    def parse_data_line(self, line : Line) -> None:
        columns = line.as_columns()

        hapfile_type_str : str | None = columns.get(0)
        hapfile_type                  = self.get_associated_hapfile_type_from_str(hapfile_type_str)

        if hapfile_type == None:
            line.err()
            error_expected(f"one of {list(HapFile.CHARACTER_TYPE_ASSOCIATIONS.keys())} when defining data", hapfile_type_str, line.number)
            return

        self.store_line_into_matrix(hapfile_type, line)
        self.validate_line_in_matrix(hapfile_type, line)


    def store_line_into_matrix(self, hftp : int, line : Line):
        columns     : Columns                = line.as_columns()
        type_header : list[tuple[str, type]] = self.header[hftp]

        if (columns.count != len(type_header)):
            line.warn()

            hftpci : int = list(HapFile.CHARACTER_TYPE_ASSOCIATIONS.values()).index(hftp)
            hftpc  : str = list(HapFile.CHARACTER_TYPE_ASSOCIATIONS.keys())[hftpci]

            error_expected(f"{len(type_header)} columns for \"{hftpc}\" entry", str(columns.count), line.number)

        matrix = self.tensor[hftp]
        matrix.append(columns.content)


    def validate_line_in_matrix(self, hftp : int, line : Line):
        columns     : Columns                = line.as_columns()
        type_header : list[tuple[str, type]] = self.header[hftp]

        if line.is_flawed():
            self.skip_due_to_errors_emmited(line)
            return

        for i, col in enumerate(type_header):
            entry = columns.content[i]

            if not is_convertible_to_type(col[1], entry):
                error_expected(f"a(n) {str(col[1])[7:-1]} for column \"{col[0]}\"", entry, line.number)
                line.warn()



    #
    # Extra
    #


    def skip_due_to_errors_emmited(self, line : Line) -> None:
        print(f">>>>> Skipping line #{line.number}.", file = stderr)



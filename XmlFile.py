import xml.etree.ElementTree as et


class XmlFile:
    def __init__(self, file_path):
        self.file_path: str = file_path
        self.start: str = ""
        self.s1: list = []    # spectrum/probe 1
        self.s2: list = []    # spectrum/probe 2
        self.probe1_pattern: str = ""       # for length, may be useful for something else
        self.probe2_pattern: str = ""
        self.seq_len: int = 0       # expected final length

        self.read_file()

    def read_file(self):
        try:
            tree = et.parse(self.file_path)
        except IOError:
            print(f"Problem reading file {self.file_path}")
            raise

        root = tree.getroot()
        self.start = root.attrib["start"]
        self.seq_len = int(root.attrib["length"])

        self.probe1_pattern = root[0].attrib["pattern"]
        self.probe2_pattern = root[1].attrib["pattern"]

        for oligon in root[0]:
            self.s1.append(oligon.text)

        for oligon in root[1]:
            self.s2.append(oligon.text)

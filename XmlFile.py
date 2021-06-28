import xml.etree.ElementTree as et


class XmlFile:
    def __init__(self, file_path):
        self.file_path = file_path
        self.start = None
        self.s1 = []    # spectrum 1
        self.s2 = []    # spectrum 2
        self.read_file()

    def read_file(self):
        try:
            tree = et.parse(self.file_path)
        except IOError:
            print(f"Problem reading file {self.file_path}")
            raise

        root = tree.getroot()
        self.start = root.attrib["start"]

        for oligon in root[0]:
            self.s1.append(oligon.text)

        for oligon in root[1]:
            self.s2.append(oligon.text)

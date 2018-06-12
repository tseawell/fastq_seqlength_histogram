import argparse
import matplotlib.pyplot as plt

comment_delimiter = "//"
sequence_delimiter = "@"
qual_score_delimiter = "+"

class fastq(object):
    id = ""
    description = ""
    sequence_str = ""
    qual_str = ""

    def __init__(self, id, description, sequence_str, qual_str):
        self.id = id
        self.description = description
        self.sequence_str = sequence_str
        self.qual_str = qual_str

    def __repr__(self):
        return "%s%s %s\n%s\n%s\n%s" %(sequence_delimiter, self.id, self.description, self.sequence_str, qual_score_delimiter, self.qual_str)

    def __str__(self):
        return self.__repr__()

    def len(self):
        return len(self.sequence_str)

def fastq_generator(infastq):
    header = ""
    sequence_str = ""

    infile = open(infastq, 'r', encoding='utf-8')

    with infile:
        is_sequence = True
        is_qual = False
        for line in infile:
            line = line.strip()
            if not line: continue
            if not is_qual and line.startswith(comment_delimiter): continue

            if not is_qual and line.startswith(sequence_delimiter):
                header = line[1:]
            elif line == qual_score_delimiter:
                is_qual = True
                is_sequence = False
            else:
                if is_sequence:
                    sequence_str = line
                    is_sequence = False
                else:
                    is_qual = False
                    yield __create_fastq__(header, sequence_str, line)
                    is_sequence = True

def open_file(infastq):
    if infastq.endswith("fastq"):
        return open(infastq, 'r', encoding='utf-8')
    else:
        print("fastq_generator requires a fastq file as input")
        exit()

def __create_fastq__(header, sequence_str, qual_str):
    header_tokens = list(header.split(" ", 2))
    id = header_tokens[0]
    description = ""
    if len(header_tokens) > 1:
        description = header_tokens[1]

    return fastq(id, description, sequence_str, qual_str)

def histogram(infile):
    print("Producing read length histogram for fastq file '%s'" %infile)
    x = []
    for seq in fastq_generator(infile):
        x.append(seq.len())

    print("File contained %d sequences" %len(x))
    num_bins = max(x)

    plt.hist(x, num_bins, facecolor='blue', alpha=0.5)
    plt.yscale('log', nonposy='clip')
    plt.xlabel('Sequence length')
    plt.ylabel('Number sequences')
    plt.title("%s" %infile)
    plt.show()


def main(arg1):
    histogram(arg1)

if __name__== "__main__":
    parser = argparse.ArgumentParser(description='Produces histogram of fastq sequence lengths.')
    parser.add_argument('infile_name', help="fastq file name")
    args = parser.parse_args()
    main(args.infile_name)



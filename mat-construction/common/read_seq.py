def read_seq(file):
    fileOpen = open( file, "r" )
    data = fileOpen.readlines()
    fileOpen.close()

    seqs = dict()

    for line in data:

        if( line[0]==">"):
            lineDivision = line.split("\n")
            seq = lineDivision[0][1:]
            seqs[seq] = ""

        else:
            lineDivision = line.split("\n")
            seqs[seq] += lineDivision[0]
    
    return (seqs)

        

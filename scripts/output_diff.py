


#This is a quick script to look at the difference between 
def output_diff(f1,f2):
    OUT1 = {}
    MATCHED = []
    DIFF = []
    with open(f2, "r") as f: 
        for line in f: 
            if line in OUT1:
                sys.exit("This shouldn't happen: only unique files should be here")
            else:
                OUT1[line] = 1
    f.close()


    with open(f1, "r") as f: 
        for line in f: 
            if line in OUT1:
                MATCHED.append(line)
            else:
                DIFF.append(line)
    f.close()

    for i in DIFF:
        print(i[:-1])


if __name__ == "__main__":
    import sys 
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    output_diff(f1,f2)

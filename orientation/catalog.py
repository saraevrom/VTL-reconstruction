

MARKUP = {
    "Star": slice(None, 8),
    "Comp": slice(8, 9),
}



if __name__=="__main__":
    with open("catalog.dat", "r") as src:
        with open("catalog.csv", "w") as dst:
            dst.write(",".join(list(MARKUP.keys()))+"\n")
            for entry in src:
                fields = list(
                    map(
                        lambda x: entry[x].strip(),
                        MARKUP.values()
                    )
                )
                dst.write(",".join(fields)+"\n")

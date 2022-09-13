def filterInPlasmoDBorESTs(inputBed, inputPlasmoDBbed, inputESTbed, outputYesBed, outputNoBed):
    """Writes only those junctions in plasmoDB or ESTs out to the "outputYesBed" 
    and only those not in ESTs or PlasmoDB into the outputNoBed."""

    ests = readJunctionsFromBed(inputESTbed)
    pdb = readJunctionsFromBed(inputPlasmoDBbed)

    outYes = open(outputYesBed, "w")
    outNo = open(outputNoBed, "w")

    for line in open(inputBed):
        if line.startswith("track"):
            continue

        pieces = line.split()
        chr = pieces[0]
        leftEdge, rightEdge = getEdges(int(pieces[1]), pieces[10], pieces[11])
        isEST = hasJunction(ests, chr, leftEdge, rightEdge)
        isPDB = hasJunction(pdb, chr, leftEdge, rightEdge)

        if isEST or isPDB:
            outYes.write(line)
        else:
            outNo.write(line)


def scoreHistogram(inputBed, binSize=10):
    """Creates the data for a score histogram."""

    d = {}
    for line in open(inputBed):
        if line.startswith("track"):
            continue

        pieces = line.split()
        score = float(pieces[4])
        bin = round(score / binSize, 0)
        if not d.has_key(bin):
            d[bin] = 0
        d[bin] += 1

    allk = d.keys()
    allk.sort()
    for k in allk:
        print k*10, "=", d[k]
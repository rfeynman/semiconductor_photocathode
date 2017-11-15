
def get1DDataInRange(x_min, x_max, nx) :
    x = []
    dx = (x_max - x_min)/(float(nx))
    xc = x_min
    x.append(xc)
    for n in range(nx) :
        xc += dx
        x.append(xc)
    return x

def get1DDataFromFile(fn, xindx, yindx, numLinesToScip=0,sep=None) :
    fh = open(fn)
    x = []
    y = []
    line_cntr = 0
    for line in fh.readlines() :
      line_cntr += 1
      if (line_cntr > numLinesToScip) : 
        vl = line.split(sep)
	print("line # = %5d, x = %s, y = %s" % (line_cntr, vl[xindx], vl[yindx]))
        x.append(float(vl[xindx]))
        y.append(float(vl[yindx]))
    fh.close()
    return (x, y)

def dump1DDataToFile(x, y, fn = "xy.dat") :
    fh = open(fn, "w")
    for i in range(len(x)) :
        fh.write("%13.5E   %13.5E\n" % (x[i], y[i]))
    fh.close()
    return

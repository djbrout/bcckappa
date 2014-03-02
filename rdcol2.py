

# Takes in Filename, reads file columnwise, and returns dictionary such that:
# import rdcol
# a = rdcol.read('filename',headline,datastartline)
# a["Column_Name"] -> returns the list for that column
#
# headline and datastartline are always > 0
#
# By Dillon Brout
# dbrout@physics.upenn.edu

def read(filename,headline,startline):
    linenum = 0
    go = 0
    column_list = []
    return_cols = {}
    inf = open(filename)
    print 'found file'
    for line in inf:
        cols = line.split()
        if linenum == headline - 1:
            for col in cols:
                return_cols[col] = []
                column_list.append(col)
                go = go + 1
        if (linenum >= startline - 1):
            index = 0
            for col in cols:
                try:
                    return_cols[column_list[index]].append(float(col))
                except:
                    return_cols[column_list[index]].append(col)
                index = index + 1
        linenum = linenum + 1
    inf.close()
    return return_cols
                

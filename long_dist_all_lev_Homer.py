import numpy as np
import matplotlib.pyplot as plt
import sys, csv, optparse
import scipy.stats

def nb_contacts_plot(matr, dom_d):
    matr = np.ma.filled(matr, 0.0)
    size = matr.shape[0]
    xy = np.array([], dtype=np.float32)
    for i in range(size):
        xy = np.append(xy,i)
    for i in range(size):
        #print i,matr[i]
        yki = matr[i]
        yki[i-min(i,1): i+min(1,size-i)] = 0
        
        if yki.sum() == 0.0 or (dom_d[i][1]-dom_d[i][0])<30: continue
        #print i,  np.argmax(yki)
        plt.suptitle('%i' %i)
        plt.plot(xy,yki)
        plt.ylim(0,600)
        plt.show()
        
def domens_matr(domens, raw_matr):
    nr_dom = len(domens.keys())
    raw_matr = np.ma.masked_invalid(raw_matr)
    raw_matr = np.ma.filled(raw_matr, 0.0)
    #new_matr = np.ma.array(np.zeros((nr_dom,nr_dom), dtype=np.float))
    new_matr = np.ma.array(np.zeros((nr_dom,nr_dom)))
    indeksy = np.arange(len(domens.keys()))
    for i in domens.keys():
        for j in domens.keys():
            #print i, j, domens[i][0], domens[i][1],domens[j][0], domens[j][1]
            #print len(raw_matr[domens[i][0]:domens[i][1]+1]), len(raw_matr[domens[j][0]:domens[j][1]+1])
            #new_matr[i,j]= raw_matr[domens[i][0]:domens[i][1]+1, domens[j][0]:domens[j][1]+1].sum()/(len(raw_matr[domens[i][0]:domens[i][1]+1]) * len(raw_matr[domens[j][0]:domens[j][1]+1]))
            #print domens[i], domens[j], new_matr
            new_matr[domens[i][0],domens[j][0]]= raw_matr[domens[i][1]:domens[i][2]+1, domens[j][1]:domens[j][2]+1].sum()
            indeksy[domens[i][0]] = i
    #new_matr = np.nan_to_num(new_matr)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #print domens.keys(), raw_matr[domens[i][1]:domens[i][2]+1, domens[j][1]:domens[j][2]+1]
    bars = ax1.bar(indeksy, +np.array(map(int, new_matr[domens['1449'][0]])), color='blue', edgecolor='black')
    bars = ax1.bar(indeksy, -np.array(map(int, new_matr[domens['1474'][0]])), color='red', edgecolor='black')
    for x,y in zip(indeksy,map(int, new_matr[domens['1449'][0]])):
        if y > 3000.0:
            plt.text(x+0.4, y+0.05, '%.2f' % y, ha='center', va= 'bottom')
    
    plt.ylim(-10000.0,+10000.0)
    plt.show()
    #sys.exit(1)
    
    return new_matr

def clip(arr, stddevs=10):
  arr = np.ma.masked_invalid(arr)
  mean = np.mean(arr)
  stddev = np.var(arr) ** 0.5
  np.clip(arr, 0, mean + stddevs * stddev, out=arr)
  return arr
  
def zero_for_beans(mat_ze, bean_nr):
    l = mat_ze.shape[0]
    ro = l
    co = 0
    for di in range(l):
        np.fill_diagonal(mat_ze[0:ro, co:l], 0.0)
        ro -= 1
        co += 1
        if co == bean_nr+1: break
    mat_ze = np.triu(mat_ze).T + np.triu(mat_ze)
    return mat_ze

def calculate_prob(mat_dom, dom):
    print "#domain_1 domain_2 p-value"
    #N = np.sum(mat_dom)
    #print N0
    N =  np.sum(np.triu(mat_dom))
    #print N
    #print np.sum(np.diagonal(mat_dom))
    acc = 0
    dl = mat_dom.shape[0]
    all_cont = dl*(dl-1)/2
    #print dom.keys()
    for  i in sorted(dom.keys()):
        for j in sorted(dom.keys()):
            #print "I_J", type(i), j
            a = sum(mat_dom[dom[i][0]][:])
            b = sum(mat_dom[dom[j][0]][:])
            #p = 1- scipy.stats.hypergeom(round(N), round(a), round(b)).cdf(round(mat_dom[i][j]))
            p = scipy.stats.hypergeom(round(N), round(a), round(b)).sf(round(mat_dom[dom[i][0]][dom[j][0]]))
            oczekiw = float(a)*float(b)/float(N)
            #if i == '1449' and j == '1474': print a,b,N,round(mat_dom[dom[i][0]][dom[j][0]])
            #print "Oczekiw",oczekiw, mat_dom[dom[i][0]][dom[j][0]]
            if oczekiw == 0.0 or mat_dom[dom[i][0]][dom[j][0]] == 0.0: 
                print i, j, 1.0
                continue
#==============================================================================
            #print "oczekiwana ", oczekiw, " obserwowana ", mat_dom[i][j]
            #print i, dom[i][0], dom[i][1], j, dom[j][0], dom[j][1], N, a, b,mat_dom[i][j], p, acc, all_cont
            #print i, dom[i][1], dom[i][2], j, dom[j][1], dom[j][2], p
            print i, j, p
#==============================================================================
            #if p < 0.1:
            #    acc +=1
            #    print "oczekiwana ", oczekiw, " obserwowana ", mat_dom[i][j]
            #    print i, j, N, a, b,mat_dom[i][j], p, acc, all_cont
    #print "Odsetek = ", float(acc)/float(all_cont)
#==============================================================================

def symmetric(ma):
    new = np.zeros((ma.shape[0], ma.shape[0]))
    for c in range(ma.shape[0]):
        new[c] = ma[c]+ma[:,c]
        new[:,c] = ma[c]+ma[:,c]
    #new[np.isnan(new)] = 0
    #print  np.allclose(np.transpose(new), new), new.size - np.isnan(new).sum()
    return new


def int_wrapper(read):
    for v in read:
        yield map(int,v[2:])

if __name__ == '__main__':

    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-m', type = "string", default = "", dest="Matrix", help = "Numpy matrix in npy format")
    optparser.add_option('-d', type = "string", default = "", dest="Domains", help ="Txt file with domain information")
    optparser.add_option('-c', type = "string", default = "", dest="Chrom", help ="The level of sherpa")
    optparser.add_option('-b', type = "int", default = 0, dest="Bin", help ="The bean number for zero values")

    (opts,args) = optparser.parse_args()
    #print len(sys.argv)
    if len(sys.argv) < 7:
        print "No matrix or domain information were given, sorry. Run script by: python long_dist.py matrix.npy domain_info.txt"
        print optparser.format_help()
        sys.exit(1)


    #if len(sys.argv) < 3:
    #    print "No matrix or domain information were given, sorry. Run script by: python long_dist.py matrix.npy domain_info.txt"
     #   sys.exit(1)
      
    else: arr = np.load(opts.Matrix)
    #print sys.argv[2] 
    #print numline
    with open(opts.Domains) as csvfile:   
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        #reader1 = int_wrapper(reader)
        #print type(reader), numline
        reader.next()
        #for r in reader:
        #   print r[0].split('chr')[1]
        #dom_dict = {rows[1]:[rows[2], rows[3]] for rows in reader}
        dom = {float(rows[1])/50000.0 :[rows[3], (float(rows[2])/50000.0)-1] for rows in reader if rows[0].split('chr')[1]== opts.Chrom}  
        
        #dom_dict = {int(k):map(int,v) for k, v in dom_dict.iteritems()}
        #dom_dict = map(int, dom_dict)
    #print "DOmeny", dom
    numline = len(dom.keys())
    dom_dict = {dom[start][0]:[nr, start, dom[start][1]] for start,nr in zip(sorted(dom.keys()), range(numline)) }
    dom.keys().sort()
    csvfile.close()    
    #print dom_dict
    
    #arr = clip(arr)
    #arr = symmetric(arr)
    #sys.exit(1)
    arr = (arr+arr.T)/2 #make arr symmetry
    bins = opts.Bin
    #print bins, opts.Bin, type(bins)
    if bins != 0:
        arr = zero_for_beans(arr, bins)
    else: pass
    #print arr[0].shape, len(dom_dict.keys())
    new_arr = domens_matr(dom_dict, arr)
    #print new_arr[0], type(new_arr[0][3]) 
        
    #nb_contacts_plot(new_arr, dom_dict)
    
    calculate_prob(new_arr, dom_dict)
    
  

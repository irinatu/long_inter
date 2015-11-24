import numpy as np
import matplotlib.pyplot as plt
import sys, csv, math, pickle, optparse
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.mlab import PCA
#import seaborn as sns

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
        print i,  np.argmax(yki)
        plt.suptitle('%i' %i)
        plt.plot(xy,yki)
        plt.ylim(0,600)
        plt.show()
        
def domens_matr(domens, raw_matr):
    nr_dom = len(domens.keys())
    #raw_matr = np.ma.masked_invalid(raw_matr)
    #raw_matr = np.ma.filled(raw_matr, 0.0)
    #new_matr = np.ma.array(np.zeros((nr_dom,nr_dom), dtype=np.float))
    new_matr = np.ma.array(np.zeros((nr_dom,nr_dom)))
    for  i in range(nr_dom):
        for j in range(nr_dom):
            #print i, j, domens[i][0], domens[i][1],domens[j][0], domens[j][1]
            #print len(raw_matr[domens[i][0]:domens[i][1]+1]), len(raw_matr[domens[j][0]:domens[j][1]+1])
            #new_matr[i,j]= raw_matr[domens[i][0]:domens[i][1]+1, domens[j][0]:domens[j][1]+1].sum()/(len(raw_matr[domens[i][0]:domens[i][1]+1]) * len(raw_matr[domens[j][0]:domens[j][1]+1]))
            new_matr[i,j]= raw_matr[domens[i][0]:domens[i][1]+1, domens[j][0]:domens[j][1]+1].sum()
            #print i,j, new_matr[i,j]
    #new_matr = np.nan_to_num(new_matr)
    
    return new_matr

def clip(arr, stddevs=10):
  arr = np.ma.masked_invalid(arr)
  mean = np.mean(arr)
  stddev = np.var(arr) ** 0.5
  np.clip(arr, 0, mean + stddevs * stddev, out=arr)
  return arr

def calculate_prob(mat_dom, dom):
    N = np.sum(mat_dom)
    acc = 0
    dl = mat_dom.shape[0]
    all_cont = dl*(dl-1)/2
    for  i in range(dl):
        for j in range(dl):
            a = sum(mat_dom[i][:])
            b = sum(mat_dom[j][:])
            #p = 1- scipy.stats.hypergeom(round(N), round(a), round(b)).cdf(round(mat_dom[i][j]))
            p = scipy.stats.hypergeom(round(N), round(a), round(b)).sf(round(mat_dom[i][j]))
            if math.isnan(p): continue
            oczekiw = float(a)*float(b)/float(N)
            if oczekiw == 0.0 or mat_dom[i][j] == 0.0: continue
#==============================================================================
            #print "oczekiwana ", oczekiw, " obserwowana ", mat_dom[i][j]
            #print i, dom[i][0], dom[i][1], j, dom[j][0], dom[j][1], N, a, b,mat_dom[i][j], p, acc, all_cont
            print i, dom[i][0], dom[i][1], j, dom[j][0], dom[j][1], p
#==============================================================================
            #if p < 0.1:
            #    acc +=1
            #    print "oczekiwana ", oczekiw, " obserwowana ", mat_dom[i][j]
            #    print i, j, N, a, b,mat_dom[i][j], p, acc, all_cont
    #print "Odsetek = ", float(acc)/float(all_cont)
#==============================================================================

def symmetric(ma):
    ma = np.ma.masked_invalid(ma)
    ma = np.ma.filled(ma, 0.0)
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
        
def dist_normalization(ma):
    l = ma.shape[0]
    row = l
    col = 0
    print ma
    for di in range(l):
        mean = ma.diagonal(di).sum()/ma.diagonal(di).shape[0]
        if mean == 0.0:
            #print "Zero", mean, di
            continue
        #print mean, di
        #print ma.diagonal(di)
        np.fill_diagonal(ma[0:row, col:l], ma.diagonal(di)/float(mean))
        row -= 1
        col += 1
        if row == 1: break
    ma = np.triu(ma).T + np.triu(ma)
    np.fill_diagonal(ma, ma.diagonal()/2.0)
    return ma
    
def plot(array):
    plt.imshow(array,origin='lower',norm=LogNorm(), interpolation='nearest')
    plt.colorbar()
    plt.show()
    
def permutation(mac, domains, macierze):
    if len(macierze) == 0:
    #macierze = []
        l = mac.shape[0]
        mat_nr = 1000
        print "The number of arrays: ", mat_nr
        for i in range(mat_nr):
            ro = l
            co = 0
            macierze.append(np.copy(mac))
            for di in range(l):
                np.fill_diagonal(macierze[i][0:ro, co:l], np.random.permutation(macierze[i].diagonal(di)))
                ro -= 1
                co += 1
                if ro == 1: break
        name = opts.Matrix.split('.')[0] + opts.Domains.split('.')[0]+".pick"
        pickle.dump(macierze, open(name, 'w'))
    else: 
        mat_nr = len(macierze)
        print mat_nr
    nr_dom = len(domains.keys())
    for  i in range(nr_dom):
        for j in range(nr_dom): 
            medians = [np.median(m[domains[i][0]:domains[i][1]+1, domains[j][0]:domains[j][1]+1]) for m in macierze] # lista median dla oddzialywania 2 domen i i j
            our_med = np.median(mac[domains[i][0]:domains[i][1]+1, domains[j][0]:domains[j][1]+1])
            greater = [g for g in medians if g >= our_med]
            #print medians, our_med, greater
            p = float(len(greater))/float(mat_nr)
            print i, domains[i][0], domains[i][1], j, domains[j][0], domains[j][1], p
            #if p != 0.0 and p != 1.0 : print 'HURAAAAAA', p
            
def principle_component(mac):
    zmiennosc = np.std(mac, axis=0)
    if len(zmiennosc) == mac.shape[0]: 
        print zmiennosc, np.where(zmiennosc==0.0), mac.shape, type(np.where(zmiennosc==0.0)[0].tolist())
        mac1 = np.delete(mac, np.where(zmiennosc==0.0)[0].tolist(), axis =1)
        print np.std(mac1, axis=0), mac1.shape

    else: print len(zmiennosc), arr_nor.shape[0]
    
    results = PCA(mac1)
    print "results", results.fracs[0]
    print "skladowe", results.Y[0:,0], type(results.Y[0:,0])
    plt.figure(1)
    #plt.subplot(211)
    plt.subplot2grid((6,6), (0, 1), rowspan=1, colspan=3)
    plt.plot(results.Y[:,0])
    print results.Y[:,0].shape, np.min(results.Y[:,0]),np.max(results.Y[:,0])
    plt.axis([0, results.Y[:,0].shape[0], np.min(results.Y[:,0]),np.max(results.Y[:,0])])
    #plt.subplot(212)
    plt.subplot2grid((6,6), (2, 0), rowspan=5, colspan=5)
    plt.imshow(mac,origin='lower',norm=LogNorm(), interpolation='nearest')
    #plt.colorbar()
    #plt.tight_layout()
    plt.show()
    
    
    

if __name__ == '__main__':

    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-m', type = "string", default = "", dest="Matrix", help = "Numpy matrix in npy format")
    optparser.add_option('-d', type = "string", default = "", dest="Domains", help ="Txt file with domain information")
    optparser.add_option('-l', type = "string", default = "", dest="Loadmtx", help ="The pickle file with list of permutated matrices")

    (opts,args) = optparser.parse_args()
    print len(sys.argv)
    if len(sys.argv) < 4:
        print "No matrix or domain information were given, sorry. Run script by: python long_dist.py matrix.npy domain_info.txt"
        print optparser.format_help()
        sys.exit(1)
      
    else: arr = np.load(opts.Matrix)
    #print sys.argv[2] 
    if opts.Loadmtx != '': LOADMATR = True
    else: LOADMATR = False
    with open(opts.Domains) as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        #reader1 = int_wrapper(reader)
        #print reader.next()
        dom_dict = {rows[1]:[rows[2], rows[3]] for rows in reader}
        dom_dict = {int(k):map(int,v) for k, v in dom_dict.iteritems()}
        #dom_dict = map(int, dom_dict)
        
    #print dom_dict
    
    arr = clip(arr)
    #plot(arr)
    arr = symmetric(arr)
    #print np.isnan(arr)
    arr_nor = dist_normalization(arr)
    plot(arr_nor)
    
    
    #print np.isnan(arr_nor)
    #arr_nor[np.isnan(arr_nor)] = 0.0
    #print np.isnan(arr_nor)
    principle_component(arr_nor)

    #if LOADMATR:
     #   macie = pickle.load(open(opts.Loadmtx))
    #else: macie = []
    #permutation(arr_nor, dom_dict, macie)
    
    
    
    #fig = plt.figure()
    
    
   
    
    #calculate_prob(new_arr, dom_dict)
    
  

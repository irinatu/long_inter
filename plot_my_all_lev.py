from sys import argv, exit
from collections import defaultdict

import optparse, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

from util import clip_and_blur, flip_to_diagonal

from matplotlib.backends import backend_pdf
from matplotlib import rcParams, cm



def lines(path,header=True):
    with open(path,'r') as handle:
        if header:
            handle.next()
        else: pass
        for line in handle:
            yield line.split(' ')

def parse_domains(gen):
    domains = defaultdict(list)
    for l in gen:
        #print l
        start = int(l[2]) 
        end = int(l[3].split("\n")[0])
        lev = l[0] 
        domains[("lev"+lev, l[1])].append((start,end))
    return domains, int(lev)

def plot_all(mtx, doms, inp, out, n_lev):
    with backend_pdf.PdfPages("%s-%s.pdf" % (os.path.basename(inp), out.split("/")[-1])) as pdf:
        #plt.figure(dpi=1200)
        #mtx = np.load(mat)
        #fig = plt.figure(dpi=2200)
        fig = plt.figure()
        colormap = plt.cm.gist_ncar
        colors = [colormap(i) for i in np.linspace(0, 0.6, n_lev)]
        #plt.imshow(mtx,interpolation='nearest',origin='lower',norm=LogNorm(),cmap=cm.jet)
        #mtx = flip_to_diagonal(mtx)
        #plt.imshow(mtx,origin='lower',norm=LogNorm(),cmap=cm.jet, interpolation='nearest')
        #plt.colorbar()
        st = 0
        for le in range(1, n_lev+1):
            doms_for_lev = dict(filter(lambda((a,b),c): a == "lev"+str(le), doms.items()))
            xyki = []
            #print doms
            #for dom_n in doms.keys():
            col_sub = "23%i" %le
            plt.subplot(int(col_sub))
            plt.imshow(mtx,origin='lower',norm=LogNorm(),cmap=cm.jet, interpolation='nearest')
            for dom_n in range(st, st + len(doms_for_lev)):
                #print dom_n, doms[dom_n], [[(start,start),(end,start),(end,end)] for start,end in doms[dom_n]]
                #print sum([[(start,start), (end,start), (end,end)] for start, end in doms[dom_n]],[])
                print "lev"+str(le)
                xyki = xyki + sum([[(start,start), (end,start), (end,end)] for start, end in doms_for_lev[("lev"+str(le), str(dom_n))]],[])
            #xyki = xyki + sum([[(end,start), (end, end)] for start, end in doms[dom_n]],[])
###########################################################################################
            #if dom_n !=0 and dom_n%10 == 0:
            #    print dom_n, doms[str(dom_n)], doms[str(dom_n)][0][1] - doms[str(dom_n)][0][0]
            #    start = doms[str(dom_n)][0][0]
            #    end = doms[str(dom_n)][0][1]
            #    plt.text(start, start-50, dom_n, fontsize=5)
############################################################################################
                start = doms_for_lev[("lev"+str(le), str(dom_n))][0][0]
                end = doms[("lev"+str(le), str(dom_n))][0][1]
                print doms_for_lev[("lev"+str(le), str(dom_n))], start, end, dom_n, le
                if end - start > 30: # names of domains are written if domains are bigger than 30
                    #print dom_n, doms[str(dom_n)], doms[str(dom_n)][0][1] - doms[str(dom_n)][0][0]
                    plt.text(start, start-30, dom_n, fontsize=2, color= colors[le-1])
            st += len(doms_for_lev)    
            
                 
            #print dom_n, xyki
            x,y = zip(*xyki)
            #col_sub = "23%i" %le
            #plt.subplot(int(col_sub))
            #plt.plot(x,y,linestyle='-',linewidth=0.5,alpha=0.7, color = colors[le-1], label = "lev"+str(le))
            #plt.imshow(mtx,origin='lower',norm=LogNorm(),cmap=cm.jet, interpolation='nearest')
            plt.plot(x,y,linestyle='-',linewidth=0.3,alpha=0.7, color = colors[le-1])
        
        #plt.legend(bbox_to_anchor=(-0.5, 1.), loc = 'upper center')
            plt.axis([0,mtx.shape[0],0,mtx.shape[0]])
        #plt.legend(loc='best')
            plt.title("Plot %i" %le,fontsize=7)
        pdf.savefig(fig)
        plt.close()
        print "finished plotting" 

if __name__=="__main__":
    
    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-m', type = "string", default = "", dest="Matrix", help = "Numpy matrix in npy format")
    optparser.add_option('-d', type = "string", default = "", dest="Domains", help ="Txt file with domain information")
    
    (opts,args) = optparser.parse_args()
    if len(argv) ==1:
        print optparser.format_help()
        exit(1)

    domeny, num_l = parse_domains(lines(opts.Domains, header=False))
    matr = np.load(opts.Matrix)
    
    matr = clip_and_blur(matr)
    plot_all(matr, domeny, opts.Matrix, opts.Domains, num_l)
    
    
    

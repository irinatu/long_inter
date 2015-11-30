from sys import argv, exit
from collections import defaultdict

import optparse, os, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

from util import clip_and_blur

from matplotlib.backends import backend_pdf
from matplotlib import rcParams, cm



def lines(path,header=False):
    with open(path,'r') as handle:
        if header:
            handle.next()
        else: pass
        for line in handle:
            yield line.split('\n')[0].split(' ')

def parse_domains(gen):
    domains = defaultdict(list)
    for l in gen:
        #print l
        start = int(l[2]) 
        end = int(l[3].split("\n")[0])
        lev = l[0] 
        domains[("lev"+lev, l[1])].append((start,end))
    return domains, int(lev)

def plot_all(mtx, mt_i, inp, out):
    with backend_pdf.PdfPages("%s-%s.pdf" % (os.path.basename(inp), os.path.basename(out))) as pdf:
        #plt.figure(dpi=1200)
        #mtx = np.load(mat)
        #fig = plt.figure(dpi=2200)
        fig = plt.figure()
        #colormap = plt.cm.gist_ncar
        #plt.imshow(mtx,interpolation='nearest',origin='lower',norm=LogNorm(),cmap=cm.jet)
        
        plt.imshow(mtx,origin='lower',norm=LogNorm(),cmap="Blues", interpolation='nearest')
        #print "SPRAWDZAM", mt_i[1196:1219, 1640:1672]    
        plt.imshow(mt_i,origin='lower',norm=LogNorm(),cmap="Reds", alpha = 0.5)
        plt.colorbar()
        
        plt.axis([0,mtx.shape[0],0,mtx.shape[0]])
        plt.legend(loc='best')
        plt.title("Plot",fontsize=7)
        pdf.savefig(fig, dpi = 1500)
        plt.close()
        print "finished plotting" 

def prepar_interac_matr(inte, si, dom):
    i_m = np.zeros((si, si))
    ZERO = True
    #print dom
    for l in inte:
        #print l[6], "tak"
        if l[6] != '0.0':
            #print "nie 0.0"
            dom1 = [dom[i] for i in dom.keys() if i[1] == l[0]]
            dom2 = [dom[i] for i in dom.keys() if i[1] == l[3]]
            dom1_2 = sum(dom1 + dom2,[])
            
            if len(dom1_2) == 2:
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = round(-math.log10(float(l[6])), 5)
                
            elif len(dom1_2) == 1:
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = round(-math.log10(float(l[6])), 5)
                
            else: print "More domains!!!!", dom1_2
        elif ZERO:
            
            dom1 = [dom[i] for i in dom.keys() if i[1] == l[0]]
            dom2 = [dom[i] for i in dom.keys() if i[1] == l[3]]
            dom1_2 = sum(dom1 + dom2,[])
           
            if len(dom1_2) == 2:
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 1000
                
            elif len(dom1_2) == 1:
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 1000
                
            else: print "More domains!!!!", dom1_2
        else:
            pass
    
    return i_m

if __name__=="__main__":
    
    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-m', type = "string", default = "", dest="Matrix", help = "Numpy matrix in npy format")
    optparser.add_option('-d', type = "string", default = "", dest="Domains", help ="Txt file with domain information")
    optparser.add_option('-i', type = "string", default = "", dest="Interaction", help ="Txt file with domain-domain interactions")
    
    (opts,args) = optparser.parse_args()
    if len(argv) ==1:
        print optparser.format_help()
        exit(1)

    domeny, le = parse_domains(lines(opts.Domains, header=False))
    #print type(domeny)
    matr = np.load(opts.Matrix)
    
    matr = clip_and_blur(matr)
    size =  matr.shape[0]
    
    inter_mtx = prepar_interac_matr(lines(opts.Interaction, header=False), size, domeny)
    
    plot_all(matr, inter_mtx, opts.Matrix, opts.Interaction)
    
    #Sprawdzenie symetrycznosci
    #inter_mtx[np.isnan(inter_mtx)] = 0
    #print  np.allclose(np.transpose(inter_mtx), inter_mtx), inter_mtx.size - np.isnan(inter_mtx).sum()
    
    

from sys import argv, exit
from collections import defaultdict

import optparse, os, math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

from util import clip_and_blur

from matplotlib.backends import backend_pdf
from matplotlib import rcParams, cm



def lines(path,header=True):
    with open(path,'r') as handle:
        if header:
            handle.next()
        else: pass
        for line in handle:
            yield line.split('\n')[0].split(' ')

def parse_domains(gen):
    domains = defaultdict(list)
    end0 = 0.0
    for l in gen:
        #print l
        if l[1] == opts.Level and l[2] == opts.Chrom:
            start = float(l[3])/150000.0 
            end = (float(l[4].split("\n")[0])/150000.0) - 1
            lev = l[1]
            if end > end0: end0 = end 
            domains[("lev"+lev, l[0])].append((start,end))
    return domains, end0

def plot_all(mtx, mt_i, out):
    with backend_pdf.PdfPages("%s_roznicowa_Rafala.pdf" % ( os.path.basename(out))) as pdf:
        #plt.figure(dpi=1200)
        #mtx = np.load(mat)
        #fig = plt.figure(dpi=2200)
        fig = plt.figure()
        #colormap = plt.cm.gist_ncar
        #plt.imshow(mtx,interpolation='nearest',origin='lower',norm=LogNorm(),cmap=cm.jet)
        
        plt.imshow(np.tril(mtx),origin='lower',norm=LogNorm(),cmap="Blues", interpolation='nearest')
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
    header = True
    if header: inte.next()
    else: pass
    int_dict = {}
    for l in inte:
        #print l
        if l[2] != '0.0' and float(l[2]) < 0.05:
            #print "nie 0.0"
            dom1 = [dom[i] for i in dom.keys() if i[1] == l[0]]
            dom2 = [dom[i] for i in dom.keys() if i[1] == l[1]]
            dom1_2 = sum(dom1 + dom2,[])
            #print dom1_2
            if dom1_2[0]==dom1_2[1]:
                #print dom1_2[0], dom1_2[1]
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 0.0
            
            elif len(dom1_2) == 2:
                #print dom1_2[0][0], int(dom1_2[0][1])+1, dom1_2[1][0], int(dom1_2[1][1])+1, i_m.shape
                #i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = round(-math.log10(float(l[2])), 5)
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 1.0
                domeny = "%s %s" %(l[0], l[1])
                int_dict[domeny]=l[2]
                
            elif len(dom1_2) == 1:
                #print dom1_2, dom1, dom2
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 1.0
                domeny = "%s %s" %(l[0], l[1])
                int_dict[domeny]=l[2]
                
            else: print "More domains!!!!", dom1_2
        elif l[2] == '0.0' and ZERO:
            
            dom1 = [dom[i] for i in dom.keys() if i[1] == l[0]]
            dom2 = [dom[i] for i in dom.keys() if i[1] == l[1]]
            dom1_2 = sum(dom1 + dom2,[])
           
            if dom1_2[0]==dom1_2[1]:
                #print dom1_2[0], dom1_2[1]
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 0.0
                
            elif len(dom1_2) == 2:
                #print dom1_2, type(dom1)
                #i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 500
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 1.0
                domeny = "%s %s" %(l[0], l[1])
                int_dict[domeny]=l[2]
                
            elif len(dom1_2) == 1:
                #print dom1_2, dom1, dom2
                i_m[int(dom1_2[0][0]):int(dom1_2[0][1])+1, int(dom1_2[1][0]):int(dom1_2[1][1])+1] = 1.0
                domeny = "%s %s" %(l[0], l[1])
                int_dict[domeny]=l[2]
                
            else: print "More domains!!!!", dom1_2
        else:
            pass
    
    return np.triu(i_m), int_dict

if __name__=="__main__":
    
    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-m', type = "string", default = "", dest="Matrix", help = "Numpy matrix in npy format")
    optparser.add_option('-d', type = "string", default = "", dest="Domains", help ="Txt file with domain information")
    optparser.add_option('-i', type = "string", default = "", dest="Interaction", help ="Txt file with domain-domain interactions")
    optparser.add_option('-f', type = "string", default = "", dest="Files", help ="Txt file with information abouth the whole path of domin and interactions files in the collumns in the format domain_path interaction_path ")
    optparser.add_option('-l', type = "string", default = "", dest="Level", help ="The level of sherpa")
    optparser.add_option('-c', type = "string", default = "", dest="Chrom", help ="The level of sherpa")
    
    (opts,args) = optparser.parse_args()
    if len(argv) ==1:
        print optparser.format_help()
        exit(1)

    domeny, dl = parse_domains(lines(opts.Domains, header=False))
    #print domeny
    
    inter_mtx, interac_d = prepar_interac_matr(lines(opts.Interaction, header=False), dl, domeny)
    
    matr = np.load(opts.Matrix)
    
    matr = clip_and_blur(matr)
    
    other_mtx=lines(opts.Files, header=True)
    int_mtres = []
    for lin in other_mtx:
        dom, d = parse_domains(lines(lin[0], header=False))
        int_m, dict = prepar_interac_matr(lines(lin[1], header=False), d, dom)
        int_mtres.append(int_m)
    
    #plt.imshow(inter_mtx,origin='lower',norm=LogNorm(),cmap="Reds")
    #plt.show()
    #plt.imshow(int_mtres[0],origin='lower',norm=LogNorm(),cmap="Reds")
    #plt.show()
    
    for ma in int_mtres:
        inter_mtx[np.where(inter_mtx==ma)] = 0.0
    
    for ke in domeny.keys():
        for ka in domeny.keys():
            #print domeny[ke][0][0], domeny[ka][0][0]
            if np.sum(inter_mtx[int(domeny[ke][0][0]):int(domeny[ke][0][1])+1, int(domeny[ka][0][0]):int(domeny[ka][0][1])+1]) != 0.0:
                key = "%s %s" %(ke[1], ka[1])
                #print key, np.sum(inter_mtx[int(domeny[ke][0][0]):int(domeny[ke][0][1])+1, int(domeny[ka][0][0]):int(domeny[ka][0][1])+1])
                inter_mtx[np.where(inter_mtx[int(domeny[ke][0][0]):int(domeny[ke][0][1])+1, int(domeny[ka][0][0]):int(domeny[ka][0][1])+1] == 1.0)] = round(-math.log10(float(interac_d[key])), 5)
                print ke[1], ka[1], interac_d[key]
    
    plot_all(matr, inter_mtx, opts.Interaction)
    
    
    #Sprawdzenie symetrycznosci
    #inter_mtx[np.isnan(inter_mtx)] = 0
    #print  np.allclose(np.transpose(inter_mtx), inter_mtx), inter_mtx.size - np.isnan(inter_mtx).sum()
    
    

# for extract from file with interactions of domains in the format nr_domain start_position end_position nr_domain start_position end_position p-value the interactions only for one defined level of sherpa


from sys import argv, exit
import optparse


def lines(path,header=True):
    with open(path,'r') as handle:
        if header:
            handle.next()
        else: pass
        for line in handle:
            yield line.split(' ')

def parse_domains(gen, lev, od_name):
    od_file = open(od_name, "w")
    domains = []
    for l in gen:
        #print l
        if l[0] == lev:
            domains.append(l[1])
            od_file.write(' '.join(l))
    od_file.close()         
    return domains
    
def extract(dom, inter, o_name):
    o_file = open(o_name, "w")
    for int_l in inter:
        #print int_l[0], int_l[3], int_l
        if int_l[0] in dom and int_l[3] in dom:
            o_file.write(' '.join(int_l))
            #print "TAKKK"
        else: pass
    o_file.close()


if __name__=="__main__":
    optparser = optparse.OptionParser(usage = "%prog [<options>]")
    optparser.add_option('-i', type = "string", default = "", dest="Interactions", help = "File with information about domains interactions")
    optparser.add_option('-d', type = "string", default = "", dest="Domains", help ="Txt file with domain information")
    optparser.add_option('-l', type = "string", default = "", dest="Level", help ="The sherpa's level")
    
    (opts,args) = optparser.parse_args()
    if len(argv) ==1:
        print optparser.format_help()
        exit(1)
        
    print "The level you are chosen is", opts.Level
    out_dom = opts.Domains.split('/')[-1].split('.')[0]+"_only"+opts.Level+"level.out"
    domeny = parse_domains(lines(opts.Domains, header=False), opts.Level, out_dom)
    #print domeny
    out_name = opts.Interactions.split('/')[-1].split('.')[0]+"_only"+opts.Level+"level.out"
    extract(domeny, lines(opts.Interactions, header=False), out_name)
    
         
import sys, os

def refineBPofBIM(_bim, _out):
    
    occupied={}
    with open(_bim) as fin, open(_out, 'w') as fout:
        for l in fin:
            [chr, rsid, gd, bp, al1, al2]=l.split()
            bp=int(bp)
            while bp in occupied:
                bp+=1
            occupied[bp]="%s %s %s %d %s %s\n"%(chr, rsid, gd, bp, al1, al2)
        for key in sorted(occupied):
            fout.write(occupied[key])

        
    
    return _out


if __name__ == '__main__':


    [_bim, _out] = sys.argv[1:3]

    refineBPofBIM(_bim, _out)

"""Exact ambient/Koszul cohomology and Gamma characters at one polynomial.

P1 H1 uses homogeneous Cech Laurent monomials with both exponents negative.
The coordinate interchange acts with the Cech orientation sign in H1.
No numerical rank tolerance and no generic-rank assumption are used.
"""
import sys, json
from pathlib import Path
from itertools import product, combinations
sys.path.insert(0, r'D:\Hexagon\work\eBLLM\runs\device\venv\Lib\site-packages')
import sympy as s

K=[(-1,0,0,1),(-1,-3,2,2),(0,1,-1,0),(1,1,-1,-1),(1,1,0,-2)]
C=(1,2,5,2,1)
CHARS=list(product((1,-1), repeat=2))

def deformation_terms():
    terms={}
    for coefficient,I in zip((1,2,3,5,7,11), combinations(range(4),2)):
        J=[j for j in range(4) if j not in I]
        for bits in product((0,2),repeat=2):
            ex=[0]*4
            for i in I:ex[i]=1
            for j,b in zip(J,bits):ex[j]=b
            terms[tuple(ex)]=terms.get(tuple(ex),0)+coefficient
    terms[(1,1,1,1)]=13
    return terms

G=deformation_terms()

def ambient(k):
    if -1 in k: return None, []
    q=sum(d<0 for d in k)
    # store the y exponent; the x exponent is k_i - y_i
    basis=list(product(*(range(d+1) if d>=0 else range(d+1,0) for d in k)))
    return q,basis

def eigenbasis(k, char):
    q,bs=ambient(k)
    if q is None: return []
    seen=set(); vectors=[]; lookup={b:i for i,b in enumerate(bs)}
    for i,v in enumerate(bs):
        if i in seen:continue
        j=lookup[tuple(d-e for d,e in zip(k,v))]
        seen.update([i,j])
        if (-1)**sum(v)!=char[0]:continue
        sb=(-1)**q
        if i==j:
            if sb==char[1]:vectors.append({i:1})
        else:vectors.append({i:1,j:char[1]*sb})
    return vectors

def line(k):
    k=tuple(k); km=tuple(d-2 for d in k)
    qs,src=ambient(km); qt,tgt=ambient(k)
    lookup={v:i for i,v in enumerate(tgt)}
    M=s.zeros(len(tgt),len(src)) if qs==qt and qs is not None else None
    if M is not None:
        for j,v in enumerate(src):
            for bits in product((0,1),repeat=4):
                w=tuple(vv+2*b for vv,b in zip(v,bits))
                if w in lookup:M[lookup[w],j]+=C[sum(bits)]
            # Evaluate the pencil multiplication at theta=1 ONLY to certify
            # that its determinant polynomial is nonzero.  The actual smooth
            # complex structure is theta=pi; transcendence transfers nonvanishing.
            for ex,c in G.items():
                w=tuple(vv+b for vv,b in zip(v,ex))
                if w in lookup:M[lookup[w],j]+=c
    out={str(q):{} for q in range(4)};rank_total=0
    for char in CHARS:
        se=eigenbasis(km,char);te=eigenbasis(k,char)
        rank=0
        if M is not None and se and te:
            # coefficient at a target orbit representative is the coordinate
            # of its normalized eigenvector (representative coefficient 1).
            cols=[]
            for v in se:
                col=s.zeros(len(src),1)
                for i,c in v.items():col[i]=c
                image=M*col
                coords=s.Matrix([image[min(vv)] for vv in te])
                recon=s.zeros(len(tgt),1)
                for c,vv in zip(coords,te):
                    for i,coeff in vv.items():recon[i]+=c*coeff
                assert recon==image, ('equivariance',k,char)
                cols.append(coords)
            rank=s.Matrix.hstack(*cols).rank()
        rank_total+=rank
        counts=[0]*4
        if qt is not None:counts[qt]+=len(te)
        if qs is not None:counts[qs-1]+=len(se)
        if M is not None:
            counts[qt]-=rank;counts[qs-1]-=rank
        assert all(x>=0 for x in counts)
        for q,v in enumerate(counts):out[str(q)][str(char)]=int(v)
    dims=[sum(out[str(q)].values()) for q in range(4)]
    expected=2*sum(s.prod(k[i] for i in I) for I in combinations(range(4),3))+2*sum(k)
    assert sum((-1)**q*d for q,d in enumerate(dims))==expected
    if M is not None:assert rank_total==M.rank()
    return {'degree':k,'ambient_source': [qs,len(src)], 'ambient_target':[qt,len(tgt)],
            'multiplication_rank':int(rank_total) if M is not None else None,
            'determinant_at_theta_1':str(M.det()) if M is not None and M.rows==M.cols else None,
            'dimensions':dims,'characters':out,'index':int(expected)}

def main():
    V=[line(k) for k in K]
    W=[dict(pair=[i+1,j+1],**line(tuple(a+b for a,b in zip(K[i],K[j]))))
       for i,j in combinations(range(5),2)]
    extra=s.Matrix([list(row) for row in zip(*K)]+[[1]*5]).nullspace()
    forbidden=[]
    for a in (0,2):
        for b in range(5):
            if a==b:continue
            row=line(tuple(x-y for x,y in zip(K[a],K[b])))
            assert row['dimensions'][1]==0
            forbidden.append(dict(target=a+1,source=b+1,**row))
    result={'seed_polynomial_coefficients':C,'actual_complex_structure':'theta=pi',
      'rank_certificate_evaluation_only':'theta=1',
      'deformation_monomials':[{'y_exponents':k,'coefficient':v} for k,v in G.items()],
      'character_order':CHARS,'line_bundles':V,
      'wedge2':W,'singlet_positive_charge_vanishings':forbidden,
      'extra_U1_kernel':[[str(x) for x in v] for v in extra],
      'full_GSC_proved':False,'method':'exact rational ranks, Cech Gamma eigenvectors'}
    target=Path(__file__).with_name('spectrum-results.json')
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('V h^q:',[x['dimensions'] for x in V])
    print('wedge2 (pair,h^q,rank):',[(x['pair'],x['dimensions'],x['multiplication_rank']) for x in W])
    for name,rows in [('V',V),('wedge2V',W)]:
        for q in [1,2]:print(name,'H'+str(q),[sum(r['characters'][str(q)][str(ch)] for r in rows) for ch in CHARS])
    print('extra U1:',result['extra_U1_kernel'])

if __name__=='__main__':main()

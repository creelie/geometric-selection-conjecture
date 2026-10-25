import spectrum as m
import sympy as s
from itertools import product
from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
k=(-2,-3,2,3);km=tuple(d-2 for d in k)
qs,src=m.ambient(km);qt,tgt=m.ambient(k);lookup={v:i for i,v in enumerate(tgt)}
M=s.zeros(24,24)
for j,v in enumerate(src):
 for bits in product((0,1),repeat=4):
  w=tuple(vv+2*b for vv,b in zip(v,bits))
  if w in lookup:M[lookup[w],j]+=m.C[sum(bits)]
 for ex,c in m.G.items():
  w=tuple(vv+b for vv,b in zip(v,ex))
  if w in lookup:M[lookup[w],j]+=c
cert=[];tex=[]
for ch in m.CHARS:
 se=m.eigenbasis(km,ch);te=m.eigenbasis(k,ch);cols=[]
 for v in se:
  a=s.zeros(24,1)
  for i,c in v.items():a[i]=c
  im=M*a
  cols.append(s.Matrix([im[min(vv)] for vv in te]))
 N=s.Matrix.hstack(*cols);det=N.det()
 assert N.shape==(6,6) and det!=0
 cert.append({'character':ch,'matrix':[[int(x) for x in row] for row in N.tolist()],'determinant':int(det)})
 chname=''.join('+' if x==1 else '-' for x in ch)
 tex.append('\\[\n M_{'+chname+'}= '+s.latex(N).replace('\\left[','\\left(').replace('\\right]','\\right)')
            +',\\qquad \\det M_{'+chname+'}='+str(det)+'.\n\\]\n')
productdet=s.prod(x['determinant'] for x in cert)
assert productdet==393460125696
# Source and target character basis determinants have the same total sign;
# verify against the full lexicographic-basis matrix as well.
assert M.det()==productdet
(root/'matrix-certificates.tex').write_text('\n'.join(tex),encoding='utf8')
(root/'checks/matrix-certificates.json').write_text(json.dumps(cert,indent=2)+'\n',encoding='utf8')
print('PASS: four character blocks, determinants', [x['determinant'] for x in cert])

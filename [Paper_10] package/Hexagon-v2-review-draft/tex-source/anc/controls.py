"""Independent finite controls; not a full physical GSC certificate."""
import sys,json,math
from pathlib import Path
from fractions import Fraction as Q
from itertools import product, combinations, combinations_with_replacement
sys.path.insert(0,r'D:\Hexagon\work\eBLLM\runs\device\venv\Lib\site-packages')
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
K=[(-1,0,0,1),(-1,-3,2,2),(0,1,-1,0),(1,1,-1,-1),(1,1,0,-2)]
B=[b for b in product(range(4),repeat=5) if sum(b)==5]
assert len(B)==101
canon=lambda b:min(tuple((x+j)%5 for x in b) for j in range(5))
assert len(set(map(canon,B)))==101
ct=[0]*4
for aa,bb,cc in combinations_with_replacement(B,3):
 e=tuple(a+b+c for a,b,c in zip(aa,bb,cc))
 for r in range(4):
  if all(x%5==(3-r)%5 for x in e):ct[r]+=1
assert ct==[566,425,280,160]
mass=[0]*4;partners={}
for b in B:
 found=[]
 for c in B:
  for r in range(4):
   if all((x+y+1)%5==(3-r)%5 for x,y in zip(b,c)):found.append((c,r))
 assert len(found)==1
 c,r=found[0];partners[b]=c;mass[r]+=1
assert mass==[51,30,20,0]
assert all(partners[partners[b]]==b for b in B)
assert sum(b==c for b,c in partners.items())==1

rows=[('X5',[5],[1]*5),('X6',[6],[1,1,1,1,2]),('X8',[8],[1,1,1,1,4]),
 ('X10',[10],[1,1,1,2,5]),('X33',[3,3],[1]*6),('X42',[4,2],[1]*6),
 ('X322',[3,2,2],[1]*7),('X2222',[2,2,2,2],[1]*8),
 ('X44',[4,4],[1,1,1,1,2,2]),('X66',[6,6],[1,1,2,2,3,3]),
 ('X43',[4,3],[1,1,1,1,1,2]),('X64',[6,4],[1,1,1,2,2,3]),
 ('X62',[6,2],[1,1,1,1,1,3])]
table=[];h=s.symbols('h')
for name,d,w in rows:
 assert sum(d)==sum(w)
 degree=Q(math.prod(d),math.prod(w))
 # Formal product computed by truncated exact convolution.
 coeff=[Q(1),Q(0),Q(0),Q(0)]
 for wi in w:
  coeff=[coeff[j]+(wi*coeff[j-1] if j else 0) for j in range(4)]
 for di in d:
  coeff=[sum(coeff[j-k]*(-di)**k for k in range(j+1)) for j in range(4)]
 assert coeff[1]==0
 c2=degree*coeff[2];chi=degree*coeff[3]
 vals=[degree*k**3/6+c2*k/12 for k in [1,2,3]]
 assert all(v.denominator==1 for v in [degree,c2,chi]+vals)
 gcd=math.gcd(*[int(v) for v in vals]);need=abs(chi)/6
 if need.denominator==1:
  factors=list(s.factorint(int(need)))
  assert any(gcd%p!=0 for p in factors)
 else:factors=[]
 table.append({'type':name,'degree':int(degree),'c2_L':int(c2),'Euler':int(chi),
               'chi_kL':list(map(int,vals)),'gcd':gcd,'required_group_order':str(need)})

M=s.Matrix(K).T
assert M.rank()==3
assert M*s.Matrix([-4,1,6,-4,1])==s.zeros(4,1)
indices=[2*sum(math.prod(k[i] for i in I) for I in combinations(range(4),3))+2*sum(k) for k in K]
assert indices==[0,-8,0,0,-4]
c2pairs=[]
for l in range(4):
 c2pairs.append(-sum(k[i]*k[j]*2 for k in K for i,j in combinations([x for x in range(4) if x!=l],2)))
assert c2pairs==[24,8,20,12]
result={'Dwork_labels':101,'Dwork_distinct_characters':101,'Dwork_couplings_by_power':ct,
 'Dwork_mass_entries_by_power':mass,'Dwork_holomorphic_det_power':70,
 'free_quotient_table':table,'line_indices':indices,'c2_V_pairings':c2pairs,
 'residual_U1':[-4,1,6,-4,1],'full_GSC_proved':False}
Path(__file__).with_name('control-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
print('PASS: Dwork 101 labels, 1431 couplings, mass multiplicities; 13 Chern/index rows; bundle arithmetic.')

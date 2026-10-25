"""Exact certificates in the analytic smoothness proof for the theta=0 seed.

The paper proves how these identities cover every coordinate stratum and
every multiplicity partition; this file checks the polynomial identities.
It does not substitute a finite point search for smoothness over C.
"""
import sys, json
from pathlib import Path
from itertools import combinations
sys.path.insert(0,r'D:\Hexagon\work\eBLLM\runs\device\venv\Lib\site-packages')
import sympy as s

z=s.symbols('z1:5');t,x,y=s.symbols('t x y');c=[1,2,5,2,1]
f=sum(c[len(I)]*s.prod(z[i] for i in I) for r in range(5) for I in combinations(range(4),r))
Q=lambda x,y:5+2*(x+y)+x*y
for i,j in combinations(range(4),2):
    k,l=[a for a in range(4) if a not in (i,j)]
    assert s.expand(s.diff(f,z[i])-s.diff(f,z[j])-(z[j]-z[i])*Q(z[k],z[l]))==0
cert=[]
for k,r in [(0,3),(1,3),(0,4)]:
    q=sum(s.binomial(r,j)*c[k+j]*t**j for j in range(r+1))
    a,b,d=s.gcdex(q,s.diff(q,t),t)
    assert d==1 and s.expand(a*q+b*s.diff(q,t))==1
    cert.append({'k':k,'r':r,'q':str(q),'a':str(a),'b':str(b),'identity':'a*q+b*qprime=1'})
F22=s.expand(f.subs(dict(zip(z,[x,x,y,y]))))
P=s.diff(f,z[0]).subs(dict(zip(z,[x,x,y,y])))
xx=-(2*y+5)/(y+2)
g=s.cancel(P.subs(x,xx));h=s.cancel(F22.subs(x,xx))
gn,gd=g.as_numer_denom();hn,hd=h.as_numer_denom()
a,b,d=s.gcdex(gn,hn,y)
assert d==1 and s.expand(a*gn+b*hn)==1
assert [c[k]*c[k+2]-c[k+1]**2 for k in range(3)]==[1,-21,1]
assert 5+2*(-2-2)+(-2)**2==1
assert Q(4,4)==37
assert s.rem(2+10*x+2*x*x, x*x+4*x+5,x)==2*(x-4)
result={'seed_smoothness_identity_checks':True,'face_determinants':[1,-21,1],
 'diagonal_certificates':cert,'two_pair_certificate':{
 'g':str(gn),'g_den':str(gd),'h':str(hn),'h_den':str(hd),'a':str(a),'b':str(b)},
 'fixed_point_seed_values':[1,2,5,48,16],
 'deformed_smoothness_method':'proper rational discriminant closed in A1, nonzero at theta=0; theta=pi transcendental',
 'finite_testing_claims':False,'full_GSC_proved':False}
Path(__file__).with_name('smoothness-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
print(json.dumps(result,indent=2))

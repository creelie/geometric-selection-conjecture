"""Symbolic checks of the explicitly specified local EFT, not a string lift."""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent
a,F,m1,m2=s.symbols('a F m1 m2',positive=True)
x,y,u,v,p,q=s.symbols('x y u v p q',real=True)
variables=(x,y,u,v,p,q)
potential=F**2/(1-2*a*(x*x+y*y))+m1*m1*(u*u+v*v)/2+m2*m2*(p*p+q*q)/2
origin=dict.fromkeys(variables,0)
gradient=[s.simplify(s.diff(potential,z).subs(origin)) for z in variables]
hessian=s.hessian(potential,variables).subs(origin)
assert gradient==[0]*6
assert hessian==s.diag(4*a*F**2,4*a*F**2,m1**2,m1**2,m2**2,m2**2)
t1,t2=s.symbols('t1 t2',positive=True)
n1,n2=s.symbols('n1 n2',positive=True)
kahler=-n1*s.log(t1)-n2*s.log(t2)
no_scale=sum(s.diff(kahler,t)**2/s.diff(kahler,t,2) for t in (t1,t2))
assert s.simplify(no_scale-n1-n2)==0
data={'status':'PASS','sympy_version':s.__version__,
      'local_EFT_gradient_at_origin':list(map(str,gradient)),
      'local_EFT_Hessian_diagonal':list(map(str,hessian.diagonal())),
      'no_scale_identity':str(s.simplify(no_scale)),
      'physical_heterotic_superpotential_realised':False,'GSC_proved':False}
(ROOT/'symbolic-controls.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(data,indent=2))

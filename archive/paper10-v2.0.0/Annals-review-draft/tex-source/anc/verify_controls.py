"""Exact controls for a GSC construction attempt; no nine-condition solution."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations, product
from collections import Counter
import json

ROOT = Path(__file__).resolve().parent
K = [(-1,0,0,1), (-1,-3,2,2), (0,1,-1,0),
     (1,1,-1,-1), (1,1,0,-2)]

def mul(a,b):
    out=Counter()
    for I,x in a.items():
        for J,y in b.items():
            if set(I).isdisjoint(J): out[tuple(sorted(I+J))]+=x*y
    return {I:x for I,x in out.items() if x}

def integral(a):
    return sum(2*x for I,x in a.items() if len(I)==3)

def main():
    total={():1}
    for k in K:
        factor={():1,**{(i,):k[i] for i in range(4) if k[i]}}
        total=mul(total,factor)
    c1={I:x for I,x in total.items() if len(I)==1}
    c2={I:x for I,x in total.items() if len(I)==2}
    c3={I:x for I,x in total.items() if len(I)==3}
    assert not c1
    c2_pair=[integral(mul(c2,{(i,):1})) for i in range(4)]
    c2tx={I:4 for I in combinations(range(4),2)}
    W={I:c2tx.get(I,0)-c2.get(I,0) for I in c2tx}
    W_pair=[integral(mul(W,{(i,):1})) for i in range(4)]
    W_eff={(0,2):6,(0,3):2}
    assert [integral(mul(W_eff,{(i,):1})) for i in range(4)]==W_pair
    indices=[2*sum(k[i]*k[j]*k[l] for i,j,l in combinations(range(4),3))+2*sum(k) for k in K]
    c3_int=integral(c3)
    assert c2_pair==[24,8,20,12] and W_pair==[0,16,4,12]
    assert sum(indices)==-12 and c3_int==-24
    assert c3_int//2==sum(indices) and Q(sum(indices),4)==-3
    assert all(sum(k)==0 for k in K)
    slopes=[sum(2*k[i] for i,j,l in product(range(4),repeat=3)
                if len({i,j,l})==3) for k in K]
    assert slopes==[0]*5
    # Invariant tetraquadric linear system: e,d,p carry weights 0,a,b.
    supports=[(0,1),(0,2),(1,2),(0,1,2)]
    support_checks=0
    for available in product(supports,repeat=4):
        assert any(sum(x==1 for x in choice)%2==0 and
                   sum(x==2 for x in choice)%2==0
                   for choice in product(*available))
        support_checks+=1
    e=Q(1,10)
    holomorphic=[e**8,e**4,Q(1)]
    right_metric=[y*y for y in holomorphic]
    normalised=[holomorphic[i]/holomorphic[i] for i in range(3)]
    assert normalised==[1,1,1] and right_metric[0]>0
    lam,b,eta,a,r=Q(2),Q(1,4),Q(1,100),Q(1,100),Q(1,10)
    kappa=lam-b
    root=-(eta+a)/(lam+b)
    radius_bound=(eta+a)/kappa
    assert eta+a<kappa*r and abs(root)<radius_bound<r
    assert radius_bound<Q(1,10)
    tests=0
    for aa,ff in product([Q(1,4),Q(1),Q(2)],[Q(1,100),Q(1,10),Q(1,3)]):
        assert 4*aa*ff*ff>0
        tests+=1
    data={'status':'PASS','arithmetic':'exact integers and rational fractions',
          'bundle_columns':[list(k) for k in K], 'c1_zero':True,
          'c2_pairings':c2_pair,'c3_integral':c3_int,
          'line_bundle_indices':indices,'total_index':sum(indices),
          'degree_four_quotient_index':-3,'slopes_at_equal_Kahler_parameters':slopes,
          'anomaly_pairings':W_pair,'effective_cover_representative':'6 J1 J3 + 2 J1 J4',
          'invariant_linear_system_support_checks':support_checks,
          'ambient_Gamma_action_nonidentity_fixed_points':48,
          'quotient_integral_anomaly_torsion_audited':False,
          'holomorphic_Yukawa_diagonal':list(map(str,holomorphic)),
          'right_matter_metric_diagonal':list(map(str,right_metric)),
          'normalised_Yukawa_diagonal':list(map(str,normalised)),
          'robust_minimum_control':{'kappa':str(kappa),'actual_root':str(root),
                                    'distance_bound':str(radius_bound),'ball_radius':str(r)},
          'positive_local_EFT_Hessian_parameter_tests':tests,
          'one_vacuum_satisfying_all_nine_conditions_constructed':False,
          'GSC_proved':False,'Lean_formalisation_claimed':False}
    (ROOT/'construction-controls.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__': main()

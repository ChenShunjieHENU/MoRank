"""Alpha-inverted ranking evidence from the authors' unchanged MMDAgg resamples.

The official dictionary has per-kernel p-value indexing ambiguity; do not use it.
The test decision API and all kernel/resampling operations remain unchanged.
"""
import os,sys
from pathlib import Path
for n in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[n]='1'
R=Path(__file__).resolve().parent
import numpy as np
source=(R/'vendor/mmdagg_np.py').read_text(encoding='utf-8')
original={};exec(compile(source,str(R/'vendor/mmdagg_np.py'),'exec'),original)
marker='    # Step 2: compute u_alpha_hat using the bisection method'
assert source.count(marker)==1
instrumented=source.replace(marker,'    return MMD_original, M1_sorted, M2, weights, len(kernel_bandwidths_l_list)\n'+marker)
capture={};exec(compile(instrumented,'mmdagg_resample_capture','exec'),capture)

def inverted_level(x,y,seed,B1=2000,B2=2000):
    observed,m1,m2,w,k=capture['mmdagg'](np.asarray(x).reshape(-1,1),np.asarray(y).reshape(-1,1),seed=int(seed),B1=B1,B2=B2)
    weights=np.tile(w,k)
    p=np.mean(m1>=observed[:,None],axis=1)
    u=float(np.min(p/weights))
    if u>=float(np.min(1/weights))-1e-12:
        return 1.0
    # Right limit at the observed p/weight threshold handles the strict quantile rule.
    u=min(np.min(1/weights),u+1e-9)
    indices=np.ceil((B1+1)*(1-u*weights)).astype(int)-1
    indices=np.clip(indices,0,B1)
    quantiles=m1[np.arange(len(weights)),indices]
    level=float(np.mean(np.max(m2-quantiles[:,None],axis=0)>0))
    return level

def smoke():
    import json,time
    rng=np.random.RandomState(20261106);rows=[];start=time.time()
    for j in range(12):
        x=rng.normal(size=15);y=rng.normal(j/8,size=15)
        if j==0:x[:]=0;y[:]=0
        if j==1:x=np.round(x);y=np.round(y)
        if j==2:y=y[:11]
        p=inverted_level(x,y,seed=100+j,B1=199,B2=199)
        for alpha in [.01,.05,.1,.25,.5]:
            expected=original['mmdagg'](x[:,None],y[:,None],alpha=alpha,seed=100+j,B1=199,B2=199)
            actual=int(p<=alpha)
            assert actual==expected,(j,alpha,p,expected)
            rows.append({'case':j,'alpha':alpha,'inverted_level':p,'official_decision':expected,'match':True})
    elapsed=time.time()-start
    print('MMDAgg official decision/inversion selfcheck',len(rows),'matches;',round(elapsed,2),'seconds')
    x=rng.normal(size=50);y=rng.normal(.3,size=50);start=time.time();p=inverted_level(x,y,seed=42)
    print('Default2000 onegene seconds',time.time()-start,'level',p,flush=True)

if __name__=='__main__':smoke()

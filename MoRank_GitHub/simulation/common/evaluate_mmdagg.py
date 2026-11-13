import sys,csv,gzip,time,os
from pathlib import Path
from mmdagg_adapter import np,inverted_level
R = Path(os.environ["MORANK_COMMON_DIR"])
worker=int(sys.argv[1]);assert 1<=worker<=24
for ident in range(1,25):
 if ident!=worker:continue
 inp=R/'data/neighbor'/f'expression_{ident:02d}.csv.gz'
 assert inp.is_file(), inp
 out=(R / "results")/f'mmdagg_neighbor_{ident:02d}.csv'
 existing=list(csv.DictReader(out.open(encoding='utf-8'))) if out.exists() else []
 with gzip.open(inp,'rt') as f:x=np.loadtxt(f,delimiter=',',skiprows=1)
 with (R/'data/neighbor'/f'scores_{ident:02d}.csv').open(encoding='utf-8-sig') as f:meta=list(csv.DictReader(f))
 start=time.time()
 with out.open('a' if existing else 'w',encoding='utf-8',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=['dataset','row_index','gene_index','truth','inverted_level'])
  if not existing:writer.writeheader()
  for i in range(len(existing),400):
   p=inverted_level(x[i,50:],x[i,:50],seed=902611060+ident*1000+i+1)
   assert np.isfinite(p) and 0<=p<=1
   writer.writerow({'dataset':ident,'row_index':i+1,'gene_index':meta[i]['gene_index'],'truth':meta[i]['truth'],'inverted_level':p})
   if (i+1)%25==0:
    f.flush();print('MMDAgg worker',worker,'dataset',ident,'gene',i+1,'elapsed',round(time.time()-start),flush=True)
print('MMDAGG_WORKER_COMPLETE',worker,flush=True)

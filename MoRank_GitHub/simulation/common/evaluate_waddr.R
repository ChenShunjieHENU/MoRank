R <- Sys.getenv('MORANK_COMMON_DIR')
Sys.setenv(WSKR_WADDR_SP_ONLY='true')
suppressPackageStartupMessages({library(data.table);library(waddR)});setDTthreads(1)
worker<-as.integer(commandArgs(TRUE)[1]);stopifnot(worker %in% 1:24)
for(id in 1:24){
 if(id!=worker)next
 input<-file.path(R,'data/neighbor',sprintf('expression_%02d.csv.gz',id));stopifnot(file.exists(input))
 meta<-fread(file.path(R,'data/neighbor',sprintf('scores_%02d.csv',id)));x<-as.matrix(fread(input))
 output<-file.path(R,'results',sprintf('waddr_neighbor_%02d.csv',id))
 existing<-if(file.exists(output))fread(output) else NULL;done<-if(is.null(existing))0L else nrow(existing);out<-if(done)split(existing,seq_len(done)) else list()
 start<-Sys.time()
 if(done<400)for(i in seq.int(done+1L,400L)){
  set.seed(802611060L+id*1000L+i)
  res<-tryCatch(wasserstein.test(x[i,51:100],x[i,1:50],method='SP',permnum=10000),error=function(e)c(pval=NA_real_))
  out[[i]]<-data.table(dataset=id,row_index=i,gene_index=meta$gene_index[i],truth=meta$truth[i],pval=unname(res['pval']),distance=unname(res['d.wass']),N_exc=unname(res['N.exc']))
  if(i%%25==0){fwrite(rbindlist(out,fill=TRUE),output);cat('waddR worker',worker,'dataset',id,'gene',i,'elapsed',round(as.numeric(difftime(Sys.time(),start,units='secs'))),'\n');flush.console()}
 }
 fwrite(rbindlist(out,fill=TRUE),output)
}
cat('WADDR_WORKER_COMPLETE',worker,'\n')

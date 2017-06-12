#chmake() builds a matrix with phage genome MCL clusters as rows and hosts as columns.
#the cluster input is the path to the dump file output from MCL. The hostmat is a 2-column matrix with phage accessions in the first column and known host in the second.

chmake<-function(cluster,hostmat)
{
	hostlist=sort(unique(hostmat[,2]))		#create list of the unique hosts in hostmat, sorted alphabetically
	pmcl=readLines(cluster)					#read in the MCL dump file
	chmat=array(0,dim=c(length(pmcl),length(hostlist)))		#set up empty cluster-host matrix
	colnames(chmat)=hostlist				#set up column names
	rownames(chmat)=paste("mcl_",c(1:length(pmcl)),sep="")	#set up row names
	for(i in 1:length(pmcl))				#for every MCL cluster...
	{
		temp=strsplit(pmcl[[i]],'\t')[[1]]	#split out the \t in the dump file
		hmatch=match(temp,hostmat[,1])		#match the genomes to their positions in the host matrix
		hosts=hostmat[hmatch,2]				#return the hosts corresponding to these genomes
		htab=table(hosts)					#summarize hosts as a table
		chmat[i,names(htab)]=as.numeric(htab)	#put the host values from the table into the cluster-host matrix
	}
	
	return(chmat)
}

#chmakeg() is similar to chmake() but with phage gene MCL clusters as rows. It additionally takes the genome-gene presence-absence matrix as input.

chmakeg<-function(cluster,hostmat,presmat)
{
	hostlist=sort(unique(hostmat[,2]))		#create list of the unique hosts in hostmat, sorted alphabetically
	gmcl=strsplit(readLines(cluster),"\t")	#read in and remove \t from the MCL dump file
	#gmcl=cluster							#option if already read in the MCL dump file
	nclust=length(gmcl)						#number of clusters
	nhost=length(hostlist)					#number of unique hosts
	chmatg=array(0,dim=c(nclust,nhost))		#set up cluster-host matrix
	rownames(chmatg)=paste("mcl_",c(1:nclust),sep="") #set up row names
	colnames(chmatg)=hostlist				#set up column names

	for(i in 1:nclust)						#for every MCL cluster...
	{
		gmatch=match(gmcl[[i]],colnames(presmat))	#where are the genes in cluster i in the pres/abs matrix?
		if(length(gmatch)>1)
		{
			psum=apply(presmat[,gmatch],1,sum)	#sum up the number of genes from cluster i are in each of the phages
			ppos=which(psum>0)					#which phage are represented?
		}
		if(length(gmatch)==1)					#what to do if only one match
		{
			ppos=which(presmat[,gmatch]==1)
		}
		pmatch=rownames(presmat)[ppos]			#name of those phage
		hmatch=match(pmatch,hostmat[,1])		#where do those phage show up in the host matrix
		hosts=hostmat[hmatch,2]					#what are their hosts?
		htab=table(hosts)						#tabulate host representation (Note: does not account for multiple genes per phage within a cluster)
		chmatg[i,names(htab)]=as.numeric(htab)	#return values to cluster-host matrix
	}

	return(chmatg)
}


#clust2pres() returns the genome-gene presence/absence matrix associated with a set of gene clusters produced by usearch. It assumes that
#the headers in the FAA file are all in the format "Accession_i" where the Accession is the genome accession and i refers to the ith gene in that genome
#Note: clust2pres also drops any 0-sum rows.

clust2pres<-function(clustdir,genomelist,trim=T)
{
	require(Matrix)
	setwd(clustdir)			#moves to directory containing the usearch clusters
	clustlist=list.files()	#vector of all files in directory
	
	if(trim==T)				#by default, clust2pres removes any gene clusters that have 2 or fewer members.
	{
		clustsize=numeric()

		for(i in 1:length(clustlist))
		{
			clustsize[i]=length(grep(">",readLines(clustlist[i])))	#counts number of sequences in each cluster
		}

		clustlist=clustlist[clustsize>2]	#only retains those with at least 3 members by default
	}

	pcmat=Matrix(0,nrow=length(genomelist),ncol=length(clustlist))	#sets up sparse matrix to be the genome-gene presence/absence matrix
	rownames(pcmat)=genomelist		
	colnames(pcmat)=clustlist

	for(i in 1:length(clustlist))			#for every usearch cluster...
	{
		temp=readLines(clustlist[i])						#read in the cluster FASTA
		accpos=grep(">",temp)								#find sequence headers
		accs=unique(unlist(strsplit(temp[accpos],">")))[-1]	#strip away the ">"s from each header
		temp2=strsplit(accs,"_")							#split by underscores
		gens=character()									#set up character vector to contain genome IDs in the cluster file
		for(j in 1:length(temp2))
		{
			gens[j]=paste(temp2[[j]][1],temp2[[j]][2],sep="_")		#return just the genome part of the header
		}

		gens=unique(gens)									#only keeps the unique set. Does not treat gene duplications in a special way
		pcmat[match(gens,genomelist),i]=1					#inserts 1s in the appropriate positions in the presence-absence matrix
	}

	droprow=which(apply(pcmat,1,sum)==0)					#finds any zero-sum rows
	if(length(droprow)>0){pcmat=pcmat[-droprow,]}
	return(pcmat)
}

#findICCC() returns the intracluster clustering coefficient for a set of MCL clusters and the associated network's adjacency matrix

findICCC<-function(clusters,adj)
{
	require(igraph)
	ccvec=nclust=numeric()		#ccvec will contain the clustering coefficient for each MCL cluster; nclust is a vector of the number of genes or genomes in a cluster 
	
	for(i in 1:length(clusters))	#for each cluster...
	{
		tempnet=graph_from_adjacency_matrix(adj[clusters[[i]],clusters[[i]]])	#creates a subgraph of just the nodes within a given MCL cluster
		ccvec[i]=transitivity(tempnet)		#calculates the clustering coefficient of that subgraph
		nclust[i]=length(clusters[[i]])		#number of genes or genomes in the subgraph
	}

	ccvec[is.na(ccvec)]=0					#sets any NA clusters to 0
	iccc=mean(ccvec[nclust>=3])				#finds the mean clustering coefficient, constrained to clusters with at least 3 members

	results=NULL

	results$ccvec=ccvec
	results$nclust=nclust
	results$iccc=iccc

	return(results)
}

#gbk2faa() returns an FAA file given a genbank flat file as input

gbk2faa<-function(filename)
{
	gbk=readLines(filename)		#reads in gbk file
	faa=character()				#will store the final FAA
	seqstarts=grep('/translation',gbk)	#finds where each sequence starts within the gbk
	for(i in 1:length(seqstarts))	#for each gene...
	{
		stopper=0
		tempseq=gbk[seqstarts[i]]	#initiates sequence
		counter=1
		while(stopper==0)			#continues adding to the sequence until encounters the next "\", indicating the end of the gene in the gbk file
		{
			tempseq=c(tempseq,gbk[seqstarts[i]+counter])
			stopper=stopper+length(grep("\"",gbk[seqstarts[i]+counter]))
			counter=counter+1
		}
		faa=c(faa,paste(">",strsplit(filename,"\\.")[[1]][1],"_",i,sep=""))		#adds header for sequence
		tempseq=paste(tempseq,collapse="")										#collapses the lines of sequence into one string
		tempseq=paste(strsplit(strsplit(tempseq,"\"")[[1]][2],"                     ")[[1]],collapse="")	#removes the '\' from the sequence
		faa=c(faa,paste(tempseq,collapse=""))									#adds the sequence to the FAA file
	}

	return(faa)
}

#hostpred() generates a vector of how often each host is associated with genes contained in a given genome. The most frequent host can be recovered
#from the results after running as names(which.max(y$host)), where y represents the variable containing the results of running hostpred().

hostpred<-function(acc,presmat,mimat,mclmat,phosts,corr=T)
{
	genelist=colnames(mclmat)
	phageid=which(rownames(presmat)==acc)				#which row in pres is the phage found
	phagegenes=colnames(presmat)[presmat[phageid,]==1]	#which genes are in that phage
	phagemclgenes=match(phagegenes,genelist)			#where are these found in the columns of the mcl matrix
	phagemcl=character()
	for(i in 1:length(phagemclgenes)){phagemcl[i]=rownames(mclmat)[mclmat[,phagemclgenes[i]]==1]}	#what are the names of the mcl clusters with the genes

	phagemclsub=match(phagemcl,rownames(mclmat))		#where do those clusters show up in the mclmatrix rows (as positions)
	phagemclsub=phagemcl[which(! is.na(phagemclsub))]	#get rid of any that are NA (can happen because of mimax trimming)
	
	hostreps=character()
	
	temprow=match(phagemclsub,rownames(mimat))	#which rows are these in the mimat
	temprow=temprow[! is.na(temprow)]			#again, get rid of any NA
	nmcl=length(temprow)			#how many mcls are left (will be used later to correct for self-info with host pred
	temphost=mimat[temprow,]					#which hosts are included
	if(length(temprow)==1){hostsum=temphost}	#if only one, record it
	if(length(temprow)>1){hostsum=apply(temphost,2,sum)}	#if more than one, report the sum
	#hostsum=apply(temphost,2,sum)
	hostreps=as.table(hostsum)					#save as table format
	
	results=NULL
	results$host=hostreps
	results$nmcl=nmcl
	if(corr==T){results=hostselfcorr(acc,results,phosts)}	#if genome already contributed to the network, this step subtracts its contribution to the host score vector.
	return(results)
}

#hostselfcorr() subtracts the contribution of a phage to its host prediction when already included in the network
	
hostselfcorr<-function(acc,hostpred,phosts)
{
	phostpos=which(phosts[,1]==acc)			#finds the row for a given accession in the hostmat
	corrhost=phosts[phostpos,2]				#returns the correspond host
	hostpredpos=which(names(hostpred$host)==corrhost)	#finds where the host occurs in the result vector within a run of hostpred() (or after if correction was off)
	newhostpred=hostpred$host				#initiates a copy of the original result vector
	newhostpred[hostpredpos]=newhostpred[hostpredpos]-hostpred$nmcl	#subtracts the number of genes in the genome (contained within the clusters being used), since each will add 1 to the known host within the result vector
	return(newhostpred)
}

#matclean() removes any all-zero rows or columns

matclean<-function(chmat)
{
	csum=apply(chmat,1,sum)			#sum up each row
	drops=which(csum==0)			#find any rows with zero sum (no hosts present)
	if(length(drops)>0)				#if there are zero-value rows, drop them
	{
		chmat=chmat[-drops,]
		csum=csum[-drops]
	}
	hsum=apply(chmat,2,sum)			#do the same going by columns (remove any hosts that no longer have representation)
	drops=which(hsum==0)			#find any zero sum columns
	if(length(drops)>0)
	{
		chmat=chmat[,-drops]		
		hsum=hsum[-drops]
	}

	return(chmat)
}

#micalc() returns the mutual information between the row and column vectors of a matrix (Here cluster-host associations). The option fullres indicates
#if micalc() should just return the value (F; default) or if it should also return the matrix of individual entries used to calcluate the value (T).

micalc<-function(chmat,fullres=F)
{
	csum=apply(chmat,1,sum)		#sum across each row (clusters)
	hsum=apply(chmat,2,sum)		#sum across each column (hosts)
	tot=sum(chmat)				#total sum of the entire matrix
	freq=chmat/tot				#represent matrix as relative to the total, rather than as absolute values
	mimat=t(t(freq/csum)/hsum)*tot*tot		#calculates the joint distribution relative to expectation if the rows and columns were independent
	mimat[mimat!=0]=log(mimat[mimat!=0])	#takes the natural log of all non-zero values. 0s are left as 0s
	mimat=freq*mimat						#multiplies by the frequency matrix
	mival=sum(mimat)						#sums all values to complete the expected value calculation

	if(fullres==T)				#return the final matrix and final value if true
	{
		results=NULL
		results$mimat=mimat
		results$mival=mival
	}

	if(fullres==F)				#return just the final value if false (default)
	{
		results=mival
	}

	return(results)
}

#mimax() randomly removes MCL clusters from the gene network and retains deletions if they improve the mutual information score without losing any genomes

mimax<-function(chmat,origpres,gmcl,nsteps)
{
	genelist=unique(unlist(gmcl))		#pulls out each unique gene (using usearch results) from the MCL dump file
	mivec=numeric()						#empty vector will contain the micalc() value after each iteration
	chmat=matclean(chmat)				#remove any 0-rows and 0-columns from the input matrix
	mivec[1]=micalc(chmat)				#calculate the first mutual information value
	pres=origpres						#initiate a presence-absence matrix as the original matrix.

	for(i in 2:nsteps)					#for the remaining steps
	{
		rnames=rownames(chmat)			#record rownames remaining in the cluster-host matrix
		tempdrop=sample(length(rnames),1)	#choose a row to drop at random
		tempmcldrop=as.numeric(strsplit(rnames[tempdrop],"mcl_")[[1]][2])	#extract the index of the associated cluster
		genematch=match(gmcl[[tempmcldrop]],genelist)	#return the position of the genes contained in that cluster within the genelist
		tempgenelist=genelist[-genematch]	#temporarily remove these genes from the genelist
		temppres=pres[,match(tempgenelist,colnames(pres))]	#temporarily drop these genes from the working presence-absence matrix
		tempsum=apply(temppres,1,sum)		#find all row sums (necessary to see if any genomes would be lost)
		if(length(which(tempsum==0))==0)	#if no genomes are lost (i.e. all rowsums greater than 0)
		{
			tempmat=chmat[-tempdrop,]		#make a temporary cluster-host matrix that drops the cluster
			tempmat=matclean(tempmat)		#clean out any new 0-columns
			tempval=micalc(tempmat)			#calculate new mutual information
			if(tempval>mivec[i-1])			#if mutual information increased, update variables to match the working copies
			{
				chmat=tempmat
				genelist=tempgenelist
				pres=temppres
			}
		}
		mivec[i]=micalc(chmat)				#save current mutual information (even if not increased)
											#Note: if any genomes are lost, then none of the temporary changes would be preserved during the next iteration
											
		if(i==nsteps/5){print("20% complete")}	#status update messages
		if(i==2*nsteps/5){print("40% complete")}
		if(i==3*nsteps/5){print("60% complete")}
		if(i==4*nsteps/5){print("80% complete")}
		if(i==nsteps){print("100% complete")}
	}

	results=NULL

	results$mat=chmat
	results$val=mivec
	results$pres=pres
	results$genelist=genelist

	return(results)
}
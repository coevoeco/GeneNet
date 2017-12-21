#hostassign() reads in a gene cluster and returns the most common host associated with the genomes harboring the genes in the cluster.

hostassign<-function(cluster,hostmat)		
{
	seqs=readLines(cluster)		#reads the cluster
	accs=seqs[grep(">",seqs)]	#finds the accessions
	accs=sapply(accs,getacc)	#truncates accessions just to the source genome
	hostmatch=hostmat[match(accs,hostmat[,1]),2]	#finds where these accessions match up in the host association table
	if(length(hostmatch)==1){hostmax=hostmatch}	#if only 1, return it
	if(length(hostmatch)>1){hostmax=names(which.max(table(hostmatch)))}	#if more than one, return the max

	return(hostmax)
}

#hostfreq() is the same as hostassign() except it returns the fraction of the genes in a cluster that are associated with the most frequent host	
hostfreq<-function(cluster,hostmat)
{
	seqs=readLines(cluster)
	accs=seqs[grep(">",seqs)]
	accs=sapply(accs,getacc)
	hostmatch=hostmat[match(accs,hostmat[,1]),2]
	maxhost=hostassign(cluster,hostmat)
	hostfreqval=length(which(hostmatch==maxhost))/length(hostmatch)

	return(hostfreqval)
}

#getacc() assumes that gene cluster FASTA files have headers in the format ">GenomeAcc_GeneNum" where GenomeAcc is a RefSeq GenBank accession for the phage genome containing genes with GeneNum being the order in which the gene of interest occurs in the particular phage genome. This function can be modified easily to account for different header formats.
getacc<-function(acc)
{
	acc=strsplit(acc,">")[[1]][2]
	acc=strsplit(acc,"_")[[1]]
	genomeacc=paste(acc[1:2],collapse="_")
	return(genomeacc)
}

#hostvectreturn() is nearly identical to hostassign() and hostfreq() except it returns the entire host table generated in hostassign(). (Future versions may simplify this format to one function with an option for which results to return, i.e. maximum result, frequency, and/or table.)
hostvecreturn<-function(cluster,hostmat)
{
	seqs=readLines(cluster)
	accs=seqs[grep(">",seqs)]
	accs=sapply(accs,getacc)
	hostmatch=hostmat[match(accs,hostmat[,1]),2]
	hostvec=table(hostmatch)
	return(hostvec)
}

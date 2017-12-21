#alignpaste() takes two alignment FASTA files with the same taxa and concatenates them into one file using the same accessions. The original files DO NOT have to be
#sorted in the same original order. They DO have to have the exact same set of accessions.

alignpaste<-function(file1,file2)
{
	
	#file1=readLines(align1)		#read in the first file in filelist
	#file2=readLines(align2)		#commented out since better is to give these in already read with a wrapper to iterate through more than 2.
	accpos1=grep(">",file1)		#determine the line numbers where the accessions are found
	accpos2=grep(">",file2)
	seqpos1=cbind(accpos1,append(accpos1[-1]-1,length(file1)))	#finds starts and stops for sequences in file1 (keeping accessions)
	seqpos2=cbind(accpos2+1,append(accpos2[-1]-1,length(file2)))	#finds starts and stop for sequences in file2 (but drops accession line)
	nseq=length(accpos1)			#count number of accessions
	accs1=character()			#prep a variable to accept the accessions
	accs2=character()
	allalign=character()
	resort=numeric()
	for(i in 1:nseq)			#iterate through number of accessions
	{
		accs1[i]=strsplit(file1[accpos1],">")[[i]][2]		#store each accession in accs1. This sets the default order for the others.
		accs2[i]=strsplit(file2[accpos2],">")[[i]][2]
	}

	for(i in 1:nseq)
	{
		resort[i]=which(accs2==accs1[i])		#finds out how the accessions are reordered in the current file relative to file1
	}

	for(i in 1:nseq)
	{
		allalign=append(allalign,append(file1[seqpos1[i,1]:seqpos1[i,2]],file2[seqpos2[resort[i],1]:seqpos2[resort[i],2]]))
	}
				
	return(allalign)
}

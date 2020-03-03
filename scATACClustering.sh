
##########################################
#Step 0: Read input & configuration files#
##########################################


if [ ${#} != 5 ] ; then
	echo -e "\nUsage: scATACClustering.sh bamfiles chromsizes configuration outdir scriptpath\n"
	echo -e "bamfiles:\tA tab delimited file where column 1 is a sample name or id column 2 is the path of the bam file and column 3 is the path of the cell barcodes.\n"
	echo -e "chromsizes:\tA tab delimited file where column 1 is the chromosome and column 2 size of the chromosome.\n"
	echo -e "configuration:\tA tab delimited file with parameters and values. Use scATACCONFIG_Default.txt as a template.\n"
	echo -e "outdir:\tThe output directory.\n"
	echo -e "scriptpath:\tThe path of this script.\n"
	exit 1
fi


BAMFILELIST=$1
CHROMSIZES=$2
CONFIGURATIONFILE=$3
OUTDIR=$4
SCRIPTPATH=$5


while IFS=$'\t' read OPTION VALUE DESCRIPTION
do
	case ${OPTION} in 
		START)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				START=${VALUE}
			fi
			;;
		END)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				END=${VALUE}
			fi
			;;
		BINSIZE)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				BINSIZE=${VALUE}
			fi
			;;
		TOPBINS)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				TOPBINS=${VALUE}
			fi
			;;
		CLUSTERRES)
			if [[ ${VALUE} =~ ^0.[[:digit:]]+$ ]]; then
				CLUSTERRES=${VALUE}
			fi
			;;
		P1MINCLUSTERSIZE)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				P1MINCLUSTERSIZE=${VALUE}
			fi
			;;
		P1NUMPCA)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				P1NUMPCA=${VALUE}
			fi
			;;
		PEAKEXT)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				PEAKEXT=${VALUE}
			fi
			;;
		TOPPEAKS)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				TOPPEAKS=${VALUE}
			fi
			;;
		P2NUMPCA)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				P2NUMPCA=${VALUE}
			fi
			;;
		TOPVARIABLEPEAKS)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				TOPVARIABLEPEAKS=${VALUE}
			fi
			;;
		MAXJOBS)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				MAXJOBS=${VALUE}
			fi
			;;
		JOBCHECKRATE)
			if [[ ${VALUE} =~ ^[[:digit:]]+$ ]]; then
				JOBCHECKRATE=${VALUE}
			fi
			;;
	esac
done < $CONFIGURATIONFILE 

echo "Configuration Summary:"
echo "START ${START}"
echo "END ${END}"
echo "BINSIZE ${BINSIZE}"
echo "TOPBINS ${TOPBINS}"
echo "CLUSTERRES ${CLUSTERRES}"
echo "P1MINCLUSTERSIZE ${P1MINCLUSTERSIZE}"
echo "P1NUMPCA ${P1NUMPCA}"
echo "PEAKEXT ${PEAKEXT}"
echo "TOPPEAKS ${TOPPEAKS}"
echo "P2NUMPCA ${P2NUMPCA}"
echo "TOPVARIABLEPEAKS ${TOPVARIABLEPEAKS}"
echo "MAXJOBS ${MAXJOBS}"
echo "JOBCHECKRATE ${JOBCHECKRATE}"


if [ -z ${START} ] || [ -z ${END} ] || [ -z ${BINSIZE} ] || [ -z ${TOPBINS} ] || [ -z ${CLUSTERRES} ] || [ -z ${P1MINCLUSTERSIZE} ] || [ -z ${P1NUMPCA} ] || [ -z ${P2NUMPCA} ] || [ -z ${PEAKEXT} ] || [ -z ${TOPPEAKS} ] || [ -z ${TOPVARIABLEPEAKS} ] || [ -z ${MAXJOBS} ] || [ -z ${JOBCHECKRATE} ];  then
	echo "Please check that the configuration variables are the correct type (int (or float for CLUSTERRES)) and complete."
	exit 1
fi


i=0
while IFS=$'\t' read -a line
do
	SAMPLEIDS[i]=${line[0]}
	BAMFILES[i]=${line[1]}
	CELLBARCODES[i]=${line[2]}
	((i++))
done < $BAMFILELIST


mkdir -p $OUTDIR
cd $OUTDIR

##############################################
#Define function for limiting child processes#
##############################################

function checkConcurrentJobs {
	#Based on the solution provided by user BruceH
	#https://stackoverflow.com/questions/6593531/running-a-limited-number-of-child-processes-in-parallel-in-bash
	while [ `jobs | wc -l` -ge $MAXJOBS ]
	do
		sleep $JOBCHECKRATE
	done
}

############################
#Step1: Generate Bin Counts#
############################
STEP=1

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 1: Generate Bin Counts" 

	STEP1OUT=${OUTDIR}/STEP1_OUTPUT

	for i in ${!BAMFILES[@]}; do
		SAMPLEDIR=${OUTDIR}/STEP1_OUTPUT/${SAMPLEIDS[i]}
		mkdir -p "${SAMPLEDIR}"
		rm -f ${SAMPLEDIR}/CompletionSummary*.txt
		#Fork multiple java processes using the trailing &
		checkConcurrentJobs; java -jar ${SCRIPTPATH}/snATACClustering.jar generatebincounts ${BAMFILES[i]} ${CELLBARCODES[i]} $CHROMSIZES $BINSIZE $SAMPLEDIR > ${STEP1OUT}/${SAMPLEIDS[i]}.out &
	done

	wait #Wait for the forked processes to finish
 
	#check that every bam file was processed correctly by looking for "CompletionSummary"
	STEPFAIL=false
	for i in ${!BAMFILES[@]}; do
		SAMPLEDIR=${OUTDIR}/STEP1_OUTPUT/${SAMPLEIDS[i]}
		if [ ! -f ${SAMPLEDIR}/CompletionSummary* ]; then
			echo "Step 1: ${SAMPLEIDS[i]} failed."
			STEPFAIL=true
		fi
	done
	if $STEPFAIL; then
		exit 1
	fi

	echo "Step 1: Complete!" 
fi


#########################
#Step2: Merge Bin Counts#
#########################
STEP=2

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 2: Merge Bin Counts"
	
	STEP2OUT=${OUTDIR}/STEP2_OUTPUT
	mkdir -p $STEP2OUT

	#create directory list
	STEP2DIRECTORYLIST=${OUTDIR}/STEP2_OUTPUT/directorylist.txt
	rm -f $STEP2DIRECTORYLIST
	touch $STEP2DIRECTORYLIST

	for i in ${!SAMPLEIDS[@]}; do
		echo "${OUTDIR}/STEP1_OUTPUT/${SAMPLEIDS[i]}" >> $STEP2DIRECTORYLIST
	done
	

	#Remove the output files if they exist
	rm -f ${STEP2OUT}/CompletionSummary*
	rm -f ${STEP2OUT}/MergedBinCounts*
	rm -f ${STEP2OUT}/TotalMergedBinCounts_*

	#Merge the bin counts

	java -jar ${SCRIPTPATH}/snATACClustering.jar mergebincounts $STEP2DIRECTORYLIST $CHROMSIZES $BINSIZE $STEP2OUT > ${STEP2OUT}/mergebincounts.out 

	#Check if completeded
	if [ ! -f ${STEP2OUT}/CompletionSummary* ]; then
		echo "Step 2 failed."
		exit 1
	fi

	echo "Step 2: Complete!"
fi


################################
#Step 3: Generate Sparse Matrix#
################################
STEP=3

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 3: Generate Sparse Matrix (Bin Counts)"

	STEP3OUT=${OUTDIR}/STEP3_OUTPUT
	mkdir -p ${STEP3OUT}
	STEP3OUTFILE=${STEP3OUT}/SparseMatrix.txt

	rm -f ${STEP3OUTFILE}
	
	#create file lists
	STEP3COUNTFILES=${STEP3OUT}/countfiles.txt
	STEP3TOTALFILES=${STEP3OUT}/totalfiles.txt

	rm -f $STEP3COUNTFILES


	ls -1 ${OUTDIR}/STEP2_OUTPUT/MergedBinCounts_* > $STEP3COUNTFILES


	rm -f $STEP3TOTALFILES
	touch $STEP3TOTALFILES
	while read line
	do
		echo "$(dirname ${line})/Total$(basename ${line})" >> $STEP3TOTALFILES
	done < $STEP3COUNTFILES

	java -jar ${SCRIPTPATH}/snATACClustering.jar generatesparsematrix $STEP3COUNTFILES $STEP3TOTALFILES $TOPBINS $STEP3OUTFILE > ${STEP3OUT}/generatesparsematrix.out
	
	#check if outfile exists
	if [ ! -f ${STEP3OUTFILE} ]; then
		echo "Step 3 failed."
		exit 1
	fi

	echo "Step 3: Complete!"
fi


###########################
#Step 4: Pass 1 Clustering#
###########################
STEP=4

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 4: Pass 1 Clustering"

	STEP3OUT=${OUTDIR}/STEP3_OUTPUT/SparseMatrix.txt

	STEP4OUT=${OUTDIR}/STEP4_OUTPUT
	mkdir -p $STEP4OUT
	

	rm -f ${STEP4OUT}/*_clusters.txt


	Rscript ${SCRIPTPATH}/snATACClustering.R ${STEP3OUT} ${P1NUMPCA} ${P1MINCLUSTERSIZE} ${CLUSTERRES} 0.8 ${STEP4OUT} > ${STEP4OUT}/Pass1Clustering.out
	

	#Check if output files exists
	STEPFAIL=false
	for i in ${!SAMPLEIDS[@]}; do
		SAMPLEFILE=${STEP4OUT}/${SAMPLEIDS[i]}_clusters.txt
		if [ ! -f $SAMPLEFILE ]; then
			echo "Step 4: ${SAMPLEIDS[i]} failed."
			STEPFAIL=true
		fi
	done
	if $STEPFAIL; then
		exit 1
	fi

	
	echo "Step 4: Complete!"
fi


#######################################
#Step 5: Split BAM files into clusters#
#######################################
STEP=5

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 5: Split Bam Files into Clusters"

	STEP4OUT=${OUTDIR}/STEP4_OUTPUT
	STEP5OUT=${OUTDIR}/STEP5_OUTPUT

	for i in ${!BAMFILES[@]}; do
		SAMPLEDIR=${STEP5OUT}/${SAMPLEIDS[i]}
		mkdir -p "${SAMPLEDIR}"
		rm -f ${SAMPLEDIR}/clusterbam_*.bam
		rm -f ${SAMPLEDIR}/CompletionSummary*.txt
		SAMPLECLUSTERFILE=${STEP4OUT}/${SAMPLEIDS[i]}_clusters.txt
		checkConcurrentJobs; java -jar ${SCRIPTPATH}/snATACClustering.jar splitbams ${BAMFILES[i]} ${CELLBARCODES[i]} $SAMPLECLUSTERFILE $SAMPLEDIR > ${SAMPLEDIR}/splitbams.out &
	done

	wait 


	STEPFAIL=false
	for i in ${!BAMFILES[@]}; do
		SAMPLEDIR=${STEP5OUT}/${SAMPLEIDS[i]}
		if [ ! -f ${SAMPLEDIR}/CompletionSummary* ]; then
			echo "Step 5: ${SAMPLEIDS[i]} failed."
			STEPFAIL=true
		fi
	done
	if $STEPFAIL; then
		exit 1
	fi

	echo "Step 5: Complete!"
fi



############################
#Step 6: Call cluster peaks#
############################
STEP=6

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 6: Peak Calling"

	
	STEP5OUT=${OUTDIR}/STEP5_OUTPUT
	STEP6OUT=${OUTDIR}/STEP6_OUTPUT

	for i in ${!SAMPLEIDS[@]}; do
		SAMPLEDIR5=${STEP5OUT}/${SAMPLEIDS[i]}
		SAMPLEDIR6=${STEP6OUT}/${SAMPLEIDS[i]}
		mkdir -p $SAMPLEDIR6

		ls -1 ${SAMPLEDIR5}/clusterbam_*.bam > ${SAMPLEDIR6}/bamfiles.txt

		MACSSUMMITS=${SAMPLEDIR6}/peakfiles.txt
		touch $MACSSUMMITS

		while read CURBAMFILE
		do
			MACSID=$(basename $CURBAMFILE | cut -d"." -f1)

			echo "${SAMPLEDIR6}/${MACSID}_summits.bed" >> $MACSSUMMITS

			checkConcurrentJobs; macs2 callpeak -t ${CURBAMFILE} -f BAMPE -n ${MACSID} -g 'hs' --nomodel --outdir "${SAMPLEDIR6}" --verbose 0 > ${SAMPLEDIR6}/${MACSID}.out &

		done < ${SAMPLEDIR6}/bamfiles.txt

	done

	wait



	STEPFAIL=false
	for i in ${!SAMPLEIDS[@]}; do
		SAMPLEDIR6=${STEP6OUT}/${SAMPLEIDS[i]}
		MACSSUMMITS=${SAMPLEDIR6}/peakfiles.txt

		while read CURPEAKFILE
		do
			if [ ! -f ${CURPEAKFILE} ]; then
				echo "Step 6: ${SAMPLEIDS[i]} failed."
				STEPFAIL=true
			fi
			
		done < ${MACSSUMMITS}
		
		wait

	done
	if $STEPFAIL; then
		exit 1
	fi


	echo "Step 6: Complete!"
fi


#######################
#Step 7: Process Peaks#
#######################
STEP=7

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 7: Process Peaks"

	STEP6OUT=${OUTDIR}/STEP6_OUTPUT
	STEP7OUT=${OUTDIR}/STEP7_OUTPUT

	mkdir -p ${STEP7OUT}

	MERGEDFILES=${STEP7OUT}/MergeBySample_files.txt
	touch $MERGEDFILES

	for i in ${!SAMPLEIDS[@]}; do
		SAMPLEDIR6=${STEP6OUT}/${SAMPLEIDS[i]}
		SAMPLEDIR7=${STEP7OUT}/${SAMPLEIDS[i]}

		mkdir -p $SAMPLEDIR7

		MACSSUMMITS=${SAMPLEDIR6}/peakfiles.txt

		echo ${SAMPLEDIR7}/MergeBySample_summits.bed >> $MERGEDFILES
		
		checkConcurrentJobs; java -jar ${SCRIPTPATH}/snATACClustering.jar processpeaks ${MACSSUMMITS} ${CHROMSIZES} ${PEAKEXT} ${SAMPLEDIR7} "MergeBySample" > ${SAMPLEDIR7}/mergebysample.out &
	
	done

	wait


	java -jar ${SCRIPTPATH}/snATACClustering.jar processpeaks ${MERGEDFILES} ${CHROMSIZES} ${PEAKEXT} ${STEP7OUT} "MergeAll" > ${STEP7OUT}/mergeall.out


	if [ ! -f ${STEP7OUT}/MergeAll_summits.bed ]; then
		echo "Step 7: Failed."
		exit 1
	fi
			

	echo "Step 7: Complete!"
fi


##############################
#Step 8: Generate Peak Counts#
##############################
STEP=8

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 8: Generate Peak Counts"

	STEP7OUT=${OUTDIR}/STEP7_OUTPUT
	STEP8OUT=${OUTDIR}/STEP8_OUTPUT
	PEAKFILE=${STEP7OUT}/MergeAll_summits.bed

	for i in ${!BAMFILES[@]}; do
		SAMPLEDIR=${STEP8OUT}/${SAMPLEIDS[i]}
		mkdir -p "${SAMPLEDIR}"
		rm -f ${SAMPLEDIR}/CompletionSummary*.txt
		checkConcurrentJobs; java -jar ${SCRIPTPATH}/snATACClustering.jar generatepeakcounts ${BAMFILES[i]} ${CELLBARCODES[i]} $PEAKFILE $PEAKEXT $SAMPLEDIR > ${SAMPLEDIR}/generatepeakscounts.out &
	done

	wait

	STEPFAIL=false
	for i in ${!BAMFILES[@]}; do
		SAMPLEDIR=${STEP8OUT}/${SAMPLEIDS[i]}
		if [ ! -f ${SAMPLEDIR}/CompletionSummary* ]; then
			echo "Step 8: ${SAMPLEIDS[i]} failed."
			STEPFAIL=true
		fi
	done
	if $STEPFAIL; then
		exit 1
	fi

	echo "Step 8: Complete!"
fi



###########################
#Step 9: Merge Peak Counts#
###########################

STEP=9

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 9: Merge Peak Counts"

	STEP8OUT=${OUTDIR}/STEP8_OUTPUT
	STEP9OUT=${OUTDIR}/STEP9_OUTPUT

	mkdir -p ${STEP9OUT}
	
	#create file lists
	STEP9COUNTFILES=${STEP9OUT}/countfiles.txt
	STEP9TOTALFILES=${STEP9OUT}/totalfiles.txt

	rm -f $STEP9COUNTFILES
	rm -f $STEP9TOTALFILES

	touch $STEP9COUNTFILES
	touch $STEP9TOTALFILES

	for i in ${!SAMPLEIDS[@]}; do
		echo "${STEP8OUT}/${SAMPLEIDS[i]}/PeakCounts_${PEAKEXT}.txt" >> $STEP9COUNTFILES
		echo "${STEP8OUT}/${SAMPLEIDS[i]}/TotalPeakCounts_${PEAKEXT}.txt" >> $STEP9TOTALFILES
	done
	

	java -jar ${SCRIPTPATH}/snATACClustering.jar mergepeakcounts $STEP9COUNTFILES $STEP9TOTALFILES $PEAKEXT $STEP9OUT > ${STEP9OUT}/mergepeakcounts.out

	
	#check if outfile exists #TODO
	if [ ! -f ${STEP9OUT}/MergedPeakCounts_${PEAKEXT}.txt ]; then
		echo "Step 9 failed."
		exit 1
	fi

	if [ ! -f ${STEP9OUT}/TotalMergedPeakCounts_${PEAKEXT}.txt ]; then
		echo "Step 9 failed."
		exit 1
	fi

	echo "Step 9: Complete!"

fi



#################################
#Step 10: Generate Sparse Matrix#
#################################
STEP=10

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 10: Generate Sparse Matrix (Peak Counts)"

	STEP9OUT=${OUTDIR}/STEP9_OUTPUT
	STEP10OUT=${OUTDIR}/STEP10_OUTPUT
	mkdir -p ${STEP10OUT}
	STEP10OUTFILE=${STEP10OUT}/SparseMatrix.txt

	rm -f ${STEP10OUTFILE}
	
	#create file lists
	STEP10COUNTFILES=${STEP10OUT}/PeakCounts.txt
	STEP10TOTALFILES=${STEP10OUT}/TotalPeakCounts.txt

	echo "${STEP9OUT}/MergedPeakCounts_${PEAKEXT}.txt" > ${STEP10COUNTFILES}
	echo "${STEP9OUT}/TotalMergedPeakCounts_${PEAKEXT}.txt" > ${STEP10TOTALFILES}


	java -jar ${SCRIPTPATH}/snATACClustering.jar generatesparsematrix $STEP10COUNTFILES $STEP10TOTALFILES $TOPPEAKS $STEP10OUTFILE > ${STEP10OUT}/step10sparsematrix.out

	
	if [ ! -f ${STEP10OUTFILE} ]; then
		echo "Step 10 failed."
		exit 1
	fi

	echo "Step 10: Complete!"

fi


############################
#Step 11: Pass 2 Clustering#
############################
STEP=11

if [ ${START} -le ${STEP} ] && [ ${END} -ge ${STEP} ]; then
	echo "Begin Step 11: Pass 2 Clustering"

	STEP10OUT=${OUTDIR}/STEP10_OUTPUT
	STEP10OUTFILE=${STEP10OUT}/SparseMatrix.txt


	STEP11OUT=${OUTDIR}/STEP11_OUTPUT
	mkdir -p $STEP11OUT
	

	rm -f ${STEP11OUT}/*_clusters.txt


	Rscript ${SCRIPTPATH}/snATACClustering.R ${STEP10OUTFILE} ${P2NUMPCA} 1 ${CLUSTERRES} 0.8 ${STEP11OUT} ${TOPVARIABLEPEAKS} > ${STEP11OUT}/Pass2Clustering.out
	

	#Check if output files exists
	STEPFAIL=false
	for i in ${!SAMPLEIDS[@]}; do
		SAMPLEFILE=${STEP11OUT}/${SAMPLEIDS[i]}_clusters.txt
		if [ ! -f $SAMPLEFILE ]; then
			echo "Step 11: ${SAMPLEIDS[i]} failed."
			STEPFAIL=true
		fi
	done
	if $STEPFAIL; then
		exit 1
	fi

	
	echo "Step 11: Complete!"
fi

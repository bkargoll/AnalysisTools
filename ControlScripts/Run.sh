#! /bin/bash
# This script run the grid jobs and then combines the output
if [ "${1}"  == "--help" ] ; then
    echo "Script to submit grid jobs (if wanted), monitor grid jobs, run Combine and purge the jobs."
    echo "Options for running this this script"
    echo " source Run.sh --Runtime <number of mintues>     Set Maximum RunTime of GRID gobs. Default 24hr. "
    echo " source Run.sh --NoCombine                       Turn off Combining the files from the GRID.  "    
    echo " source Run.sh --Submit                          Submits jobs to the GRID before starting job monitoring." 
else

    haveProxy=`voms-proxy-info --all | wc -l`;
    echo "Line in proxy: " $haveProxy
    if [ "${haveProxy}" -ge "16"  ]; then
	nmindefault=288;
	if [ "${1}"  == "--Runtime" ]; then
	    nmin=`echo "${2}/5" | bc`;
	    
	else
	    nmin=$nmindefault;
	fi
	
	if [ "${1}"  == "--Submit" ] || [ "${2}"  == "--Submit" ] || [ "${3}"  == "--Submit" ] || [ "${4}"  == "--Submit" ];   then
	    if [ -f  jobs_complete ]; then
		rm jobs_complete
	    fi
	    if [ -f jobs_submitted ]; then
		rm jobs_submitted
	    fi
	    if [ -f Set_*/out ]; then
		rm Set_*/out
	    fi
	    if [ -f Set_*/err ]; then
		rm Set_*/err
	    fi
	    
	    source Submit --SetupAndSubmit
	fi
	myruntime=`echo "${nmin}*5/60" | bc`; 
	echo "Will check jobs every 5 minutes for  ${myruntime}  hours.";
	#echo ${nmin}
	idx=0;
	while (test "$nmin" -ge "$idx" )
	  do
	  sleep 300;
	  nsets=$(ls | grep Set_ | wc -l)
	  njobs=$(cat jobs_submittedOrComplete | wc -l)
	  if [[ ${nsets} -ne ${njobs} ]]; then
	      echo "not all jobs were submitted. Retrying failed submissions..."
	      source Submit --Submit 
	  fi
	  source CheckandGet.sh  --get >& junk_CG; rm junk_CG;
	  eval=`cat jobs_submitted  | wc -l`
	  echo ${eval} " jobs still running"
	  if [[  ${eval} -eq 0 ]]; then
	      source Purge_Jobs.sh --all
	      if [ "${1}"  == "--NoCombine" ] || [ "${2}"  == "--NoCombine" ] || [ "${3}"  == "--NoCombine" ] || [ "${4}"  == "--NoCombine" ]; then
		  echo "Jobs Complete"
	      else
		  echo "Starting Combine"
		  source Combine >& log_Combine
		  echo "Job Complete"
	      fi
	      echo "finished in loop " $idx
	      let idx=nmin+1
	  else
	      let idx=idx+1 
	      echo "in loop " $idx
	  fi
	done	
	echo "Running Complete. The output of Combine has been dumped to the file log_Combine."
    else
	echo "Please setup your voms-proxy and grid certificate"
	echo "voms-proxy-init -voms cms:/cms/dcms"
	echo "grid-proxy-init"
	echo " "
	echo "For more information type: voms-proxy-info --all" 
    fi
fi

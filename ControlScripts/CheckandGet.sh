#! /bin/bash
# This script retrieves jobs from the grid and returns the numbe of jobs check and recieved
if [ "${1}"  == "--help" ]; then
    echo "Options for running this configuration file:"
    echo "--help                    Prints this message"
    echo "--detailed                Prints out the status of each job (This option goes last)"
fi
if [ "${1}"  == "--detailed" ]; then
    cat jobs_submitted | awk '{system("glite-wms-job-status " $1 )}'
fi
if [ "${2}"  == "--detailed" ]; then
    cat jobs_submitted | awk '{system("glite-wms-job-status " $1 )}'
fi


touch jobs_complete
touch jobs_get
cat jobs_submitted | awk '{system("glite-wms-job-status " $1 " | grep -i done | wc -l |awk \x27{sum++}END{if($1 == 0) print  \" glite-wms-job-output --dir "$2 " " $1 " ; mv " $2 "/$USER*/* "$2 "; rm -rf " $2 "/$USER* ; echo " $1 " " $2 "  >> jobs_complete; source " $2 "/GRIDRetrieve.sh  \"}\x27 >> jobs_get ")}'

#source jobs_get 
#rm jobs_get
echo 'Number of jobs produced:'  
cat  jobs_submitted  | wc -l
echo 'Number of jobs Complete and Downloaded:'
cat  jobs_complete  | wc -l
 

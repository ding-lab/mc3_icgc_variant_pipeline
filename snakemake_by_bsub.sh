main() {
	export PATH="/gscuser/lwang/miniconda3/envs/snake/bin:$PATH"
	submit_to_bsub 4 4000 60 ./snakemake_bsub snakemake -p all
}

echoerr() { printf "%s\n" "$*" >&2; }

submit_to_bsub() {
    # A helper function to construct neccessary constraints for a LSF job
    # For example, submit_to_bsub 4 4000 30 /path/to/log bash myscript.sh
    # means
    #   - Use 4 CPUs
    #   - Use 4000 MB memory
    #   - Runs at most 30 minutes; job will be timeout'd and killed after 30mins
    #   - Log stdout and sterr separately to /path/to/log.stdout.log and /path/to/log.stderr.log
    #   - The command `bash myscript.sh` is all the rest of the arguments.
    local n_process=$1; shift
    local memory=$1; shift
    local timeout=$1; shift
    local log_prefix=$1; shift
    local cmd=$@

    local bsub_resource="select[mem>${memory} && ncpus>=${n_process}] rusage[mem=$memory] span[hosts=1]"
    bsub -n $n_process -M ${memory}000 -W $timeout -N \
        -o ${log_prefix}.log -e ${log_prefix}.log \
        -R "'${bsub_resource}'" $cmd
}

main

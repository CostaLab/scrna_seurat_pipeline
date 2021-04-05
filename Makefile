# DO NOT RUN, modify it with your own options
submit_example_c1_c2:
	mkdir -p logs
	echo " Submitting EXAMPLE_PROJ - cond1 cond2 d1 d2 analysis"
	sbatch run_example.sh EXAMPLE_PROJ_C1
	sleep 5
	sbatch run_example.sh EXAMPLE_PROJ_C2
	sleep 5
	squeue -l | grep "vr393028"

# DO NOT RUN, modify it with your own options
submit_viz_example_c1_c2:
	echo " Submitting viz EXAMPLE_PROJ"
	sbatch run_viz_example.sh EXAMPLE_PROJ_C1
	sleep 5
	sbatch run_viz_example.sh EXAMPLE_PROJ_C2
	sleep 5
	squeue -l | grep "vr393028"

clean:
	rm logs/*_job.log

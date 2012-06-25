# @ job_name = project1
# @ comment = "simonbe@kth.se"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL;
# @ wall_clock_limit = 00:30:00
# @ notification = error
# @ notify_user = simonbe@kth.se
# @ job_type = bluegene
# @ bg_size = 8192 
# @ queue
mpirun -exe project1 -mode VN -np 32768 -verbose 1

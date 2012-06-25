# @ job_name = project1
# @ comment = "simonbe@kth.se"
# @ error = $(job_name).$(jobid).out
# @ output = $(job_name).$(jobid).out
# @ environment = COPY_ALL;
# @ wall_clock_limit = 00:59:00
# @ notification = error
# @ notify_user = simonbe@kth.se
# @ job_type = bluegene
# @ bg_size = 512
# @ queue
scalasca -analyze mpirun -exe project1 -mode VN -np 2048

import pp, numpy as np, random as rnd

def parallel_function(max_seats):
    import numpy as np, random as rnd
    seats = np.array([-1])
    for n in range(max_seats - 1):
        i = n+1
        if(i==1 or i in seats):
            d = -1
            while(d in seats):
                d = rnd.randrange(max_seats)+1
            seats = np.append(seats, d)
        else:
            seats = np.append(seats, i)
    if(float(max_seats) in seats):
        return 1
    else:
        return 0

job_server = pp.Server() 

# Define your jobs
max_k = 100000
l = []
for k in range(max_k):
    l.append(job_server.submit(parallel_function, (100,)))

# Compute and retrieve answers for the jobs.
w = 0
for k in range(max_k):
    w+=l[k]()
print(w/max_k)

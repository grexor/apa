import os
from Queue import *
from threading import *
import glob

data = glob.glob(os.path.expanduser("~/apa/data.apa/*"))

num_worker_threads = 20
q = Queue()

def worker():
    while True:
        task = q.get()
        os.system(task)
        q.task_done()

# raw
tasks = []
for lib_id in data:
    lib_id = os.path.basename(lib_id) # get lib_id from full path
    tasks.append("pybio.bamclip -bam \"~/apa/data.apa/%s/*/*/*.bam\" -image ~/apa/data.apa/%s/%s_clipping" % (lib_id, lib_id, lib_id))

for i in range(num_worker_threads):
     t = Thread(target=worker)
     t.daemon = True
     t.start()

for task in tasks:
    q.put(task)

q.join()

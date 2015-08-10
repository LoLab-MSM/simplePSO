import os
for i in range(1,101):
	os.system('python run_pso_earm.py %s >> %s.out' %(str(i),str(i)))

import pylab as plt
import numpy as np
import os


BASE = '/Users/jamespino/git/ParticleSwarmOptimization/'
os.chdir(BASE)
def extract(FOLDER):
    os.chdir(FOLDER)
    data = np.loadtxt('pso_0.txt')
    for i in range(1,250):
        if os.path.exists("pso_%s.txt" % (str(i))):
            tmp = np.loadtxt("pso_%s.txt" % (str(i)))
            data = np.column_stack((data,tmp))
    print FOLDER,np.average(data[-1,:]),np.std(data[-1,:])
    os.chdir(BASE)
    return data

directory = 'LIKE_speed_0.10_particles_10_iterations_200'
s_1_p10 = extract(directory) 
directory = 'LIKE_speed_0.10_particles_20_iterations_200'
s_1_p20 = extract(directory)  
directory = 'LIKE_speed_0.10_particles_50_iterations_200'
s_1_p50 = extract(directory)
directory = 'LIKE_speed_0.10_particles_100_iterations_200'
s_1_p100 = extract(directory)

directory = 'LIKE_speed_0.25_particles_10_iterations_200'
s_25_p10 = extract(directory)
directory = 'LIKE_speed_0.25_particles_20_iterations_200'
s_25_p20 = extract(directory)
directory = 'LIKE_speed_0.25_particles_50_iterations_200'
s_25_p50 = extract(directory)
directory = 'LIKE_speed_0.25_particles_100_iterations_200'
s_25_p100 = extract(directory)

directory = 'LIKE_speed_0.5_particles_10_iterations_200'
s_5_p10 = extract(directory)
directory = 'LIKE_speed_0.5_particles_20_iterations_200'
s_5_p20 = extract(directory)
directory = 'LIKE_speed_0.5_particles_50_iterations_200'
s_5_p50 = extract(directory)
directory = 'LIKE_speed_0.5_particles_100_iterations_200'
s_5_p100 = extract(directory)

directory = 'LIKE_speed_1.0_particles_10_iterations_200'
s_10_p10 = extract(directory)
directory = 'LIKE_speed_1.0_particles_20_iterations_200'
s_10_p20 = extract(directory)
directory = 'LIKE_speed_1.0_particles_50_iterations_200'
s_10_p50 = extract(directory)
directory = 'LIKE_speed_1.0_particles_100_iterations_200'
s_10_p100 = extract(directory)
# 
new = 'Speed0.1_particles10_iterations200'
new2 = extract(new)
tmp = np.loadtxt('pso_0.txt')
tmp1 = np.loadtxt('pso1_0.txt')
tmp2 = np.loadtxt('pso2_0.txt')
tmp3 = np.loadtxt('pso3_0.txt')
tmp4 = np.loadtxt('pso4_0.txt')
tmp5 = np.loadtxt('pso5_0.txt')
tmp7 = np.loadtxt('pso7_0.txt')
tmp8 = np.loadtxt('pso8_0.txt')
tmp9 = np.loadtxt('pso9_0.txt')
tmp10 = np.loadtxt('pso10_0.txt')
# 10 particles
#   s_1_p10, s_25_p10, s_5_p10, s_10_p10
# 20 particles
#  s_1_p20,s_25_p20,s_5_p20,s_10_p20
# 50 particles
#  s_1_p50,s_25_p50,s_5_p50,s_10_p50
# 100 particles
#  s_1_p100,s_25_p100,s_5_p100,s_10_p100
plt.loglog(s_1_p20,'g')
plt.loglog(new2,'r')
plt.loglog(tmp,'ob')
plt.loglog(tmp1,'og')
plt.loglog(tmp2,'ok')
plt.loglog(tmp3,'pink')
plt.loglog(tmp4,'red')
plt.loglog(tmp7,'blue')
plt.loglog(tmp8,'red')
plt.loglog(tmp9,'black')
plt.loglog(tmp10,'purple')
plt.show()
def plot_overlap_speed(nparticles,n10,n25,n50,n100,sub):
    
    ax = plt.subplot(2,2,sub)
    #ax= plt.subplot(111)
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    plt.plot(n10,'green',alpha=.3)
    plt.plot(n25,'red',alpha=.3)
    plt.plot(n50,'blue',alpha=.3)
    plt.plot(n100,'orange',alpha=.3) 
    plt.plot(np.average(n10,axis=1),marker='o',color='green',label='.10',linewidth=2)
    plt.plot(np.average(n25,axis=1),marker='o',color='red',label='.25',linewidth=2)
    plt.plot(np.average(n50,axis=1),marker='o',color='blue',label='.5',linewidth=2)
    plt.plot(np.average(n100,axis=1),marker='o',color='orange',label='1.',linewidth=2)
    plt.xlim(10,250)
    plt.title('Varying speed,  %s particles '%nparticles)
    plt.xlabel('Iteration number')
    plt.ylabel('Likelihood')
    plt.legend(loc=0)
    
plt.figure()    
plot_overlap_speed('10', s_1_p10, s_25_p10, s_5_p10, s_10_p10,1)
plot_overlap_speed('20',s_1_p20,s_25_p20,s_5_p20,s_10_p20,2)
plot_overlap_speed('50',s_1_p50,s_25_p50,s_5_p50,s_10_p50,3)
plot_overlap_speed('100',s_1_p100,s_25_p100,s_5_p100,s_10_p100,4)
plt.savefig('pso_speed.png',bbox_inches='tight')
plt.show()

def plot_varying_particles(n10,n20,n50,n100,speed,sub):
    #plt.figure()
    ax = plt.subplot(2,2,sub)
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    #"""
    x10 = np.linspace(0,10*200,num=200)
    x20 = np.linspace(0,20*200,num=200)
    x50 = np.linspace(0,50*200,num=200)
    x100 = np.linspace(0,100*200,num=200)
    plt.plot(x10,np.average(n10,axis=1),'green',label='10')
    plt.plot(x20,np.average(n20,axis=1),'red',label='20')
    plt.plot(x50,np.average(n50,axis=1),'blue',label='50')
    plt.plot(x100,np.average(n100,axis=1),'orange',label='100')
    """     
    plt.plot(np.average(n10,axis=1),'green',label='10')
    plt.plot(np.average(n20,axis=1),'red',label='20')
    plt.plot(np.average(n50,axis=1),'blue',label='50')
    plt.plot(np.average(n100,axis=1),'orange',label='100')
    plt.xlim(10,250)
    """
    plt.title('Varying number of particles, %s'%speed)
    plt.xlabel('Iteration number')
    plt.ylabel('Likelihood')
    plt.legend(loc=0)
    #plt.show()
plt.figure()  
plot_varying_particles(s_1_p10,s_1_p20,s_1_p50,s_1_p100,'speed=0.1',1)    
plot_varying_particles(s_25_p10,s_25_p20,s_25_p50,s_25_p100,'speed=0.25',2)
plot_varying_particles(s_5_p10,s_5_p20,s_5_p50,s_5_p100,'speed=0.5',3) 
plot_varying_particles(s_10_p10,s_10_p20,s_10_p50,s_10_p50,'speed=1.0',4) 
plt.savefig('pso_particles.png',bbox_inches='tight')
plt.show()
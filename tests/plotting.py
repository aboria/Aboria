import matplotlib.pyplot as plt
import numpy as np

for filename in ['vector_addition','daxpy','finite_difference','multiquadric']:
    data = np.loadtxt('tests/'+filename+'.csv')
    plt.figure(figsize= (6,5))
    if filename in ['vector_addition','daxpy']:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria_level1')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_level2')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
    else:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')

    plt.legend()
    plt.xlabel('N')
    plt.ylabel('MFLOPS')
    plt.savefig(filename+'.pdf')


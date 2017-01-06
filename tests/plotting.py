import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

for filename in ['vector_addition','daxpy','finite_difference','multiquadric','multiquadric_scaling']:
    data = np.loadtxt('tests/'+filename+'.csv')
    plt.figure(figsize= (6,5))
    if filename in ['vector_addition','daxpy']:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria_level1')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_level2')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
    elif filename in ['multiquadric']:
        filename_matlab = 'matlab_matrix_multiply_size'
        data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
        plt.semilogx(data_matlab[:,0],data_matlab[:,1]/1e6,label='matlab')
    elif filename in ['multiquadric_scaling']:
        filename_matlab = 'matlab_matrix_multiply_scaling'
        data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
        plt.semilogx(data_matlab[:,0],data_matlab[:,1]/1e6,label='matlab')

    else:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria_eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')

    leg = plt.legend(fancybox=True,loc='upper left')
    leg.get_frame().set_alpha(0.3)
    plt.xlabel('N')
    plt.ylabel('MFLOPS')
    plt.tight_layout()
    plt.savefig(filename+'.pdf')



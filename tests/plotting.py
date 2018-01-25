import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties
import numpy as np
import glob
import re
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)



for filename in ['benchmark_fmm','vector_addition','daxpy','finite_difference','multiquadric','linear_spring','multiquadric_scaling']:
#for filename in ['linear_spring']:
    if filename in ['benchmark_fmm']:
        files = ['tests/benchmark_fmm_1.csv','tests/benchmark_fmm_2.csv','tests/benchmark_fmm_3.csv']
        datas = []
        for f in files:
            datas.append(np.loadtxt(f))
    elif filename in ['linear_spring']:
        files = glob.glob('tests/'+filename+'*.csv')
        datas = []
        radius_div_h = []
        for f in sorted(files):
            radius_div_h.append(float(re.search('tests\/'+filename+'(.*)\.csv',f).group(1)))
            print 'reading ',f,' with r = ',radius_div_h[-1]
            datas.append(np.loadtxt(f))
    else:
        data = np.loadtxt('tests/'+filename+'.csv')

    if filename in ['benchmark_fmm']:
        f = plt.figure(figsize= (10,5))
    else:
        f = plt.figure(figsize= (6,5))

    handles = []
    if filename in ['benchmark_fmm']:
        plt.title(r'Direct summation versus FMM versus H2 (Dim=2, N=3)')
        plt.loglog(datas[1][:,0],datas[1][:,1],label='Direct')
        plt.loglog(datas[1][:,0],datas[1][:,6],label='FMM kdtree')
        plt.loglog(datas[1][:,0],datas[1][:,7],label='FMM octtree')
        plt.loglog(datas[1][:,0],datas[1][:,8],label='H2 kdtree')
        plt.loglog(datas[1][:,0],datas[1][:,9],label='H2 octtree')

        plt.loglog(datas[1][:,0],0.5e-8*datas[1][:,0]**2,linestyle='--',label='$N^2$')
        plt.loglog(datas[1][:,0],3.0e-6*datas[1][:,0],linestyle='--',label='$N$')
        #plt.loglog(datas[1][:,0],1.3e-4*datas[1][:,0],linestyle='--',label='$N$')
    elif filename in ['vector_addition']:
        plt.title(r'$a_i = b_i + c_i$')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria (Level 1)')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='Aboria (Level 2)')
        #plt.semilogx(data[:,0],data[:,3]/1e6,label='Eigen')
        plt.semilogx(data[:,0],data[:,4]/1e6,label='std::vector')
    elif filename in ['daxpy']:
        plt.title(r'$a_i = a_i + 0.1\, b_i $')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria (Level 1)')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='Aboria (Level 2)')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='Eigen')
        plt.semilogx(data[:,0],data[:,4]/1e6,label='std::vector')
    elif filename in ['multiquadric']:
        plt.title(r'$s_i = s_i+ \sum_j^N  s_j \sqrt{|\mathbf{dx}_{ij}|+c_j^2}$')
        filename_matlab = 'matlab_matrix_multiply_size'
        #data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='aboria-eigen')
        #plt.semilogx(data[:,0],data[:,3]/1e6,label='eigen')
        plt.semilogx(data[:,0],data[:,4]/1e6,label='std::vector')
        #plt.semilogx(data_matlab[:,0],data_matlab[:,1]/1e6,label='matlab')
    elif filename in ['multiquadric_scaling']:
        plt.title(r'$s_i = s_i+ \sum_j^N  s_j \sqrt{|\mathbf{dx}_{ij}|+c_j^2}$')
        filename_matlab = 'matlab_matrix_multiply_scaling'
        data_matlab = np.loadtxt('tests/'+filename_matlab+'.csv',delimiter=',')
        plt.plot(data[:,0],data[:,1]/1e6,label='Aboria')
        plt.plot(data[:,0],data[:,2]/1e6,label='Aboria-Eigen')
        plt.plot(data[:,0],data[:,3]/1e6,label='Eigen')
        plt.plot(data[:,0],data[:,4]/1e6,label='std::vector')
        plt.plot(data_matlab[:,0],data_matlab[:,1]/1e6,label='Matlab')
    elif filename in ['linear_spring']:
        plt.title(r'$dr_i = \sum_{|\mathbf{dx}_{ij}|<2r} \frac{(2r-|\mathbf{dx}_{ij}|}{|\mathbf{dx}_{ij}|}\mathbf{dx}_{ij}$')
        dashed_line = mpl.lines.Line2D([], [], color='black',linestyle='--',label='std::vector + Gromacs search')
        solid_line = mpl.lines.Line2D([], [], color='black',linestyle='-',label='Aboria')
        handles.append(dashed_line)
        handles.append(solid_line)
        for data,radius in zip(datas,radius_div_h):
            handle, = plt.loglog(data[:,0],data[:,1]/1e6,label='radius = %2.1f'%radius)
            plt.loglog(data[:,0],data[:,2]/1e6,color=handle.get_color(),linestyle='-.',label='radius = %2.1f'%radius)
            handles.append(handle)
            plt.loglog(data[:,0],data[:,3]/1e6,color=handle.get_color(),linestyle='--',label='gromacs %f'%radius)
        plt.ylim(10,1e4)
    else:
        plt.semilogx(data[:,0],data[:,1]/1e6,label='Aboria')
        plt.semilogx(data[:,0],data[:,2]/1e6,label='Aboria-Eigen')
        plt.semilogx(data[:,0],data[:,3]/1e6,label='Eigen')



    if filename in ['benchmark_fmm']:
        # Shrink current axis by 20%
        print 'asdfasdfsdafd'
        fontP = FontProperties()
        fontP.set_size('small')
        box = plt.gca().get_position()
        plt.gca().set_position([box.x0, box.y0, box.width * 0.65, box.height])
        leg = plt.legend(fancybox=True,prop=fontP,loc='center left',bbox_to_anchor=(1.0,0.5))
    elif filename in ['multiquadric_scaling']:
        leg = plt.legend(fancybox=True,loc='upper right')
    elif filename in ['multiquadric']:
        leg = plt.legend(fancybox=True,loc='lower right')
    elif filename in ['linear_spring']:
        leg = plt.legend(handles=handles,fancybox=True,loc='upper left')
    elif filename in ['daxpy','vector_addition']:
        leg = plt.legend(fancybox=True,loc='upper right')
    else:
        leg = plt.legend(fancybox=True,loc='upper left')

    leg.get_frame().set_alpha(0.3)
    plt.xlabel(r'$N$')
    if filename in ['benchmark_fmm']:
        plt.ylabel(r'$T_e$')
    elif filename in ['vector_addition', 'daxpy']:
        plt.ylabel(r'$N / T_e$ ($\times 10^6$)')
    elif filename in ['linear_spring']:
        plt.ylabel(r'$N^2 / T_e$ ($\times 10^6$)')
    else:
        plt.ylabel(r'$N^2 / T_e$ ($\times 10^6$)')
    #plt.tight_layout()
    plt.savefig(filename+'.pdf')
    plt.savefig(filename+'.svg')

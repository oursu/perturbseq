import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

#this method taken from Dixit et al., 2016, https://github.com/asncd/MIMOSCA/blob/master/GBC_CBC_pairing/fit_moi.ipynb
def moi(adata_here,pref='perturb',level='guide',gridsize=100,maxk=10):
    
    import scipy
    from numpy import unravel_index
    
    print('Computing MOI and detection probability using code from Dixit et al., 2016')

    if pref+'.'+level+'.perturbations_per_cell' not in adata_here.obs:
        print('ERROR: missing adata.obs['+pref+'.'+level+'.perturbations_per_cell], please run perturb.pp.perturbations_per_cell first')
        exit
    if pref+'.'+level+'.cells_per_perturbation' not in adata_here.uns:
        print('ERROR: missing adata.obs['+pref+'.'+level+'.cells_per_perturbation], please run perturb.pp.cells_per_perturbation first')
        exit
        
    moi_dist=np.array(list(adata_here.obs[pref+'.'+level+'.perturbations_per_cell']))
    num_virus=adata_here.uns[pref+'.'+level+'.cells_per_perturbation'].shape[0]

    n,bins=np.histogram(moi_dist,range(int(maxk)+1))

    #maximum number of viruses possible (per cell)
    #maxk

    #total number of unique barcodes
    print('number of distinct perturbations',num_virus)

    #gridsize for performing lambda and alpha search
    nums=gridsize

    #specify start and finishing MOI to search over, it is set to 0.1 and 3 here
    mois=np.linspace(0.1,5,nums) #(0.1,2,nums)
    #specify start and finishing detection probability to search over, it is set to 0.1 and 0.99 here
    detects=np.linspace(0.1,0.99,nums)

    #initialize search array
    LL=np.zeros((nums,nums))


    #loop through square grid of different poission parameters and detection probabilities
    for i in range(nums):
        for m in range(nums):

            #current parameter guesses
            moi_guess=mois[i]
            detect_guess=detects[m]

            #initialize possion distribution with current guess    
            pdf=scipy.stats.poisson.pmf(k=range(maxk),mu=moi_guess)

            #Zero truncation and renormalization
            pdf[0]=0
            pdf=np.divide(pdf,np.sum(pdf))


            #get probabilities after convolving with binomial distribution
            zibpdf=np.zeros((maxk,1))
            for k in range(maxk):
                pf=0
                for j in np.arange(k,maxk):
                    pf+=pdf[j]*scipy.stats.binom.pmf(k,j,detect_guess)
                zibpdf[k]=pf

            #evaluate log likelihood after multiplying with observed values
            ll=1.0
            for k in range(maxk):#range(len(n)):
                ll+=n[k]*np.log(zibpdf[k])
            LL[i,m]=ll

    #Log likelihood vs. paramter space
    plt.contour(np.round(detects,2),np.round(mois,2),LL,400,cmap='magma')
    plt.colorbar()
    plt.xlabel('Detection Probability')
    plt.ylabel('MOI')

    #Find parameters that maximize the log likelihood
    final_tuple=unravel_index(LL.argmax(), LL.shape)
    moi_guess=int(100*mois[final_tuple[0]])/100
    detect_guess=int(100*detects[final_tuple[1]])/100
    print('MOI:',moi_guess)
    print('Detection probability:',detect_guess)

    adata_here.uns['MOI']=moi_guess
    adata_here.uns['Detection_probability']=detect_guess
    plt.scatter(detect_guess,moi_guess,color='black',s=50)

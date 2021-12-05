#!/usr/bin/python

#imports
import numpy as np
from numpy import random as r
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition, datasets
from sklearn.cross_decomposition import CCA
from sklearn.cluster import KMeans
import scipy.stats as st
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf


#read in and subset data
data = pd.read_csv('./TEST_output.csv')




#############################################################################
# NOTE: EVERYTHING BELOW IS FROM MY OLD EPIPHYTE WORK AND NEEDS TO BE TWEAKED
#       TO WORK FOR THIS ANALYSIS
#############################################################################




#create standardized dataframe
stand = (morph - morph.mean())/morph.std()
#df of each sample average distance from means of its vars (i.e. its 'characteristic deviation')
dev = stand.mean(axis = 1)
#var means
means = morph.mean(axis=0)
#var stds
stds = morph.std(axis=0)


#write an imputation function
def impute_val(i,j):
    morph.ix[i,j] = means[j] + dev[i]*stds[j]
    return

for i,j in zip(i_list, j_list):
    impute_val(i,j)

#Now, drop rows 14 and 34, because they were missing too much
morph = morph.drop(34, axis = 0)
morph = morph.drop(14, axis = 0)




#subset the data for the environ vars of interest (and retaining only rows that had no morph NAs
env = data.ix[morph.index,['alt', 'c3m', 'c3s', 'c9m', 'c9s', 'ppt', 'slp', 'asp', 'rug', 'ndv', 'brt', 'grn', 'wet']]
#double-check no NAs in the environ data
assert env.shape == env.dropna().shape

#center and standardize data
stand_morph = (morph-morph.mean())/morph.std()
stand_env = (env-env.mean())/env.std()



#define function that returns proper axis label based on morphological column name
def get_morph_axis_label(col_name):
    if col_name.lower().startswith('dens'):
        label =  '%s %s density, %s ($\mu m^{-2}$)'
    else:
        label = '%s %s area, %s ($\mu m^{2}$)'
    if col_name[7] == 'B':
        side = 'Abaxial'
    else:
        side = 'Adaxial'
    if col_name[9] == 'B':
        region = 'basal'
    else:
        region = 'mid-leaf'
    if col_name[11] == 'T':
        struct = 'trichome'
    else:
        struct = 'stomate'
    return label % (side, struct, region)


# HISTOGRAMS
############


stand_morph_subset = stand_morph.ix[:,range(10)]
morph_subset = morph.ix[:,range(10)]
morph_subset.hist(bins = 25, xrot=45, xlabelsize=8, color = 'green')

#plt.tight_layout()

env.hist(bins = 25, xrot = 45, xlabelsize = 8, color = 'green')


##############
# SCATTERPLOTS
##############

#plot morph vars against alt
fig_morph_v_alt = plt.figure()
plt.suptitle('Morphological Vars vs. Altitude')
morph_a_d_cols = [col for col in morph.columns if (col.startswith('AREA') or col.startswith('DENS'))]
morph_a_d = morph.ix[:, morph_a_d_cols]
for i in range(morph_a_d.shape[1]):
    ax = fig_morph_v_alt.add_subplot(3,np.ceil((morph_a_d.shape[1]/3.)),i+1)
    plot(data.ix[list(morph_a_d.index), 'alt'], morph_a_d.ix[:,i], 'o', color = 'green')
    #ax.set_title(morph_a_d.columns[i])
    ax.set_xlabel('Altitude (m)')
    ax.set_ylabel(get_morph_axis_label(morph_a_d.columns[i]))
    plt.xticks(rotation=45)
    
    #NOTE: FOR SOME REASON, THIS ONLY SEEMS TO PRODUCE A NICE EFFECT IF I RUN IT IN THE COMMAND LINE WHILE A PLOT IS OPEN, RATHER THAN IF I RUN IT AS PART OF A SCRIPT...
    #plt.tight_layout() 



#plot morph vars against my exposure index
fig_morph_v_exp = plt.figure()
plt.suptitle('Morphological vals vs. Exposure')
for i in range(morph.shape[1]):
    fig_morph_v_exp.add_subplot(3,np.ceil((morph.shape[1]/3)),i+1)
    plot(data.ix[list(morph.index), 'exp'], morph.ix[:,i], 'og')
    plt.xticks(rotation=45)


#####
# PCA
#####

#Briefly playing with PCA of morpho data, following sklearn's example here:
#http://scikit-learn.org/stable/auto_examples/decomposition/plot_pca_iris.html
#pca = decomposition.PCA(n_components = 3)
#
#pca.fit(morph)
#
#tf_morph= pca.transform(morph)
#
#fig3 = plt.figure(1, figsize=(4, 3))
#plt.clf()
#ax = Axes3D(fig3, rect=[0, 0, .95, 1], elev=48, azim=134)
#plt.cla()
#
#ax.scatter(tf_morph[:, 0], tf_morph[:, 1], tf_morph[:, 2], s = 180, c=data.ix[list(morph.index), 'alt'], 
#        cmap=plt.cm.Greens, edgecolor='k')
#
#




# Alternative PCA take: using the hard-coded PCA algorithm from Laurel's lecture

#define PCA function
def PCA(A, numpc = 0):
    M = (A - mean(A.T, axis = 1)).T #subtract mean along cols
    [latent, coeff] = linalg.eig(cov(M)) #compute eigenvalues and eigenvectors of cov mat
    p = size(coeff, axis = 1) #get number of vars
    idx = argsort(latent)[::-1] #get correctly descending order of eigenvalues 
    coeff = coeff[:,idx] #sort coefficients (i.e. loadings) accordingly
    latent = latent[idx] #sort eigenvalues accordingly
    if numpc <p and numpc >= 0:
        coeff = coeff[:, range(numpc)] #take only desired num of PCs
    score = dot(coeff.T, M) #calculate scores
    return(coeff, score, latent)




#get alt and exp values, for coloring PCA points
alt_color = data.ix[list(morph.index), 'alt']
exp_color = data.ix[list(morph.index), 'exp']
exp_color = [2.5 if np.isnan(e) else e for e in exp_color]





#morph_coeff, morph_score, morph_latent = PCA(np.array(morph_subset), numpc = morph_subset.shape[1])  #dropping the binary fruit_flow and the leaf-count vars


#fig5 = plt.figure()
#
#fig5.add_subplot(231)
#plt.scatter(morph_score[0], morph_score[1], s = 120, c = alt_color, cmap = plt.cm.Greens, edgecolor = 'k')
#
#fig5.add_subplot(232)
#plt.scatter(morph_score[0], morph_score[2], s = 120, c = alt_color, cmap = plt.cm.Greens, edgecolor = 'k')
#
#fig5.add_subplot(233)
#plt.scatter(morph_score[1], morph_score[2], s = 120, c = alt_color, cmap = plt.cm.Greens, edgecolor = 'k')
#
#ax1 = fig5.add_subplot(235, projection = '3d')
#Axes3D.scatter(ax1, xs = morph_score[0], ys = morph_score[1], zs = morph_score[2], s = 120, c = alt_color, cmap = plt.cm.Greens, depthshade = False, edgecolor = 'k')
#
#ax2 = fig5.add_subplot(236, projection = '3d')
#Axes3D.scatter(ax2, xs = morph_score[0], ys = morph_score[1], zs = morph_score[2], s = 120, c = exp_color, cmap = plt.cm.Greens, depthshade = False, edgecolor = 'k')




#Plot PC-transformed data and loadings
#PCs = [0,1]
#sc = 5
#
#
#fig6 = plt.figure()
#ax = fig6.add_subplot(111)
#plt.scatter(morph_score[0], morph_score[1], s = 120, c = alt_color, cmap = plt.cm.Greens, edgecolor = 'k')
#lines = [[(0,0), (sc*coeff[0], sc*coeff[1])] for coeff in morph_coeff[:,PCs]]
#lines = mc.LineCollection(lines, colors = 'black')
#ax.add_collection(lines)
#ax.autoscale()
#[ax.text(*item) for item in zip(1.1*sc*morph_coeff[:,0], 1.1*sc*morph_coeff[:,1], cols)]



#PCAs of both datasets, for scree plots and input to CCA
m_coeff, m_score, m_latent = PCA(morph_subset, numpc = len(morph_subset.columns))
e_coeff, e_score, e_latent = PCA(stand_env, numpc = len(env.columns))

#scree plots
fig_scree = plt.figure()
fig_scree.suptitle('Scree plots of morpho and env PCA')
ax = fig_scree.add_subplot(121)
ax.plot(range(len(morph_subset.columns)), m_latent/sum(m_latent), '-o', color = 'black', linewidth = 4)
ax.xaxis.set_label_text('Morpho PCs')
ax = fig_scree.add_subplot(122)
plt.plot(range(len(stand_env.columns)), e_latent/sum(e_latent), '-o', color = 'black', linewidth = 4)
ax.xaxis.set_label_text('Env PCs')



#Plot PC-transformed data and loadings
PCs = [0,1]
sc = 5

m_cols = list(morph_subset.columns)
e_cols = list(stand_env.columns)



fig_m_loads = plt.figure()
ax = fig_m_loads.add_subplot(111)
ax.scatter(m_score[0], m_score[1], s = 120, c = e_score[0], cmap = plt.cm.Greens, edgecolor = 'k')
lines = [[(0,0), (sc*coeff[0], sc*coeff[1])] for coeff in m_coeff[:,PCs]]
lines = mc.LineCollection(lines, colors = 'black')
ax.add_collection(lines)
ax.autoscale()
[ax.text(*item) for item in zip(1.1*sc*m_coeff[:,0], 1.1*sc*m_coeff[:,1], m_cols)]
fig_m_loads.suptitle('Morpho PCs 1 and 2, with loadings, colored by env PC 1', size = 25)
ax.xaxis.set_label_text('morpho PC 1', size = 15)
ax.yaxis.set_label_text('morpho PC 2', size = 15)

fig_e_loads = plt.figure()
ax = fig_e_loads.add_subplot(111)
ax.scatter(e_score[0], e_score[1], s = 120, c = m_score[0], cmap = plt.cm.Greens, edgecolor = 'k')
lines = [[(0,0), (sc*coeff[0], sc*coeff[1])] for coeff in e_coeff[:,PCs]]
lines = mc.LineCollection(lines, colors = 'black')
ax.add_collection(lines)
ax.autoscale()
[ax.text(*item) for item in zip(1.1*sc*e_coeff[:,0], 1.1*sc*e_coeff[:,1], e_cols)]
fig_e_loads.suptitle('Env PCs 1 and 2, with loadings, colored by morpho PC 1', size = 25)
ax.xaxis.set_label_text('env PC 1', size = 15)
ax.yaxis.set_label_text('env PC 2', size = 15)




#####
# LMs
#####

#Add PCs to the two dataframes
for i in range(4):
    env['env_pc%i' % (i+1)] = e_score[i]

for i in range(4):
    morph['morph_pc%i' % (i+1)] = m_score[i]

#then combine into one df, for ease of running linear models
df = env.merge(morph, left_index = True, right_index = True)



#define function to run a regression on given y and x, get the results, and plot the fit line on currently
#active axes with solid line if significant (as long as no interaction terms)
def fit_and_plot_ols(y, x, data, ax, bonf_factor = 1):
    #run regression, plot regression line and stats
    formula = '%s ~ %s' % (y, x)
    mod = smf.ols(formula, data = data).fit()
    p = mod.pvalues[1]
    slope = mod.params[1]
    intercept = mod.params[0]
    if p < 0.05/bonf_factor: #Bonferroni correction
        linetype = '-'
    else:
        linetype = ':'
    #xmin, xmax = ax.get_xlim()
    #plot_x = np.linspace(xmin, xmax, 100)
    plt.plot(data[x], [i*slope + intercept for i in data[x]], '%sr' % linetype, linewidth = 2)
    return

pc = 'env_pc1'
#plot morph vals against 1st environ PC
fig_morph_v_env_pc1 = plt.figure()
fig_morph_v_env_pc1.suptitle('Morphological Variables vs. Environmental PC1', size = 20)
for i in range(morph_a_d.shape[1]):
    ax = fig_morph_v_env_pc1.add_subplot(3,np.ceil((morph_a_d.shape[1]/3.)),i+1)
    plot(df[pc], morph_a_d.ix[:,i], 'o', color = 'green')
    #ax.set_title(morph_a_d.columns[i])
    ax.set_xlabel('Environmental PC1')
    ax.set_ylabel(get_morph_axis_label(morph_a_d.columns[i]))
    plt.xticks(rotation=45)
    #run regression, plot regression line and stats
    fit_and_plot_ols(morph_a_d.columns[i], pc, df, ax, bonf_factor = morph_a_d.shape[1])

#plt.tight_layout()

pc = 'env_pc2'
#plot morph vals against 1st environ PC
fig_morph_v_env_pc2 = plt.figure()
fig_morph_v_env_pc2.suptitle('Morphological Variables vs. Environmental PC2', size = 20)
for i in range(morph_a_d.shape[1]):
    ax = fig_morph_v_env_pc2.add_subplot(3,np.ceil((morph_a_d.shape[1]/3.)),i+1)
    plot(df[pc], morph_a_d.ix[:,i], 'o', color = 'green')
    #ax.set_title(morph_a_d.columns[i])
    ax.set_xlabel('Environmental PC2')
    ax.set_ylabel(get_morph_axis_label(morph_a_d.columns[i]))
    plt.xticks(rotation=45)
    #run regression, plot regression line and stats
    fit_and_plot_ols(morph_a_d.columns[i], pc, df, ax, bonf_factor = morph_a_d.shape[1])

#plt.tight_layout()






#Combine all trich data into repeated-measures




#####
# CCA
#####



#------------------- USING FIRST 4 PCs FROM EACH DATASET

n_pcs = 4
n_components = 3

cca = CCA(n_components = n_components)
cca.fit(e_score[0:n_pcs].T, m_score[0:n_pcs].T)

#get transformed X and Y vars
e_cc, m_cc = cca.transform(e_score[0:n_pcs].T, m_score[0:n_pcs].T)


#text's proportional placements along x and y from lr corner (integer increments represent 0.0, 0.1, ..., 0.9)
text_x_prop = 3
text_y_prop = 1

#plot
fig7 = plt.figure()
for i in range(n_components):
    ax = fig7.add_subplot(200 + 10*n_components +1 + i)
    plt.scatter(e_cc[:,i], m_cc[:,i], s = 75, c = 'chartreuse', edgecolor = 'k')
    plt.text(np.linspace(*ax.get_xlim(),num = 10)[-text_x_prop],np.linspace(*ax.get_ylim(),num = 10)[text_y_prop], '$r^2$ = %0.3f' % np.corrcoef(e_cc[:,i].T, m_cc[:,i].T)[1,0])

#plot loadings

#create dataframe for env barplot
env_var = ['e_pc_'+str(i+1) for i in range(n_pcs)]*n_components
comp_num = np.hstack([[i]*(len(env_var)/n_components) for i in range(n_components)])
env_loads = pd.DataFrame({'env':env_var, 'loadings': cca.x_loadings_.flatten(), 'n_comp':comp_num})

#crate dataframe for morph barplot
morph_var = ['m_pc_'+str(i+1) for i in range(n_pcs)]*n_components
comp_num = np.hstack([[i]*(len(morph_var)/n_components) for i in range(n_components)])
morph_loads = pd.DataFrame({'morph':morph_var, 'loadings': cca.y_loadings_.flatten(), 'n_comp':comp_num})

#make barplots
ax = fig7.add_subplot(200 + 10*n_components + i + n_components-1)
b1 = sns.barplot(x='morph', y = 'loadings', hue = 'n_comp', data = morph_loads)
b1.set_xticklabels(b1.xaxis.get_ticklabels(), rotation=50)
ax = fig7.add_subplot(200 + 10*n_components + i + n_components + 1)
b2 = sns.barplot(x='env', y = 'loadings', hue = 'n_comp', data = env_loads)
b2.set_xticklabels(b2.xaxis.get_ticklabels(), rotation=0)



#------------------- USING ALL DATA


#num of CCs to keep
n_components = 2

#create the CCA object and fit it to the data
cca = CCA(n_components = n_components)
cca.fit(env, morph_subset)

#get transformed X and Y vars
env_c, morph_c = cca.transform(env, morph_subset)


#text's proportional placements along x and y from lr corner (integer increments represent 0.0, 0.1, ..., 0.9)
text_x_prop = 3
text_y_prop = 1

#plot
fig7 = plt.figure()
for i in range(n_components):
    ax = fig7.add_subplot(200 + 10*n_components +1 + i)
    plt.scatter(env_c[:,i], morph_c[:,i], s = 75, c = 'chartreuse', edgecolor = 'k')
    plt.text(np.linspace(*ax.get_xlim(),num = 10)[-text_x_prop],np.linspace(*ax.get_ylim(),num = 10)[text_y_prop], '$r^2$ = %0.3f' % np.corrcoef(env_c[:,i].T, morph_c[:,i].T)[1,0])



#plot loadings

#create dataframe for env barplot
env_var = list(env.columns)*n_components
comp_num = np.hstack([[i]*len(env.columns) for i in range(n_components)])
env_loads = pd.DataFrame({'env':env_var, 'loadings': cca.x_loadings_.flatten(), 'n_comp':comp_num})

#crate dataframe for morph barplot
morph_var = list(morph_subset.columns)*n_components
comp_num = np.hstack([[i]*len(morph_subset.columns) for i in range(n_components)])
morph_loads = pd.DataFrame({'morph':morph_var, 'loadings': cca.y_loadings_.flatten(), 'n_comp':comp_num})

#make barplots
ax = fig7.add_subplot(200 + 10*n_components + i + n_components)
b1 = sns.barplot(x='morph', y = 'loadings', hue = 'n_comp', data = morph_loads)
b1.set_xticklabels(b1.xaxis.get_ticklabels(), rotation=50)
ax = fig7.add_subplot(200 + 10*n_components + i + n_components + 1)
b2 = sns.barplot(x='env', y = 'loadings', hue = 'n_comp', data = env_loads)
b2.set_xticklabels(b2.xaxis.get_ticklabels(), rotation=0)




##################
# CLUSTER ANALYSIS
##################

#from stats class example I wrote on  10/10/17
def clust_plot(a,b, C_var, C_cov):
    C = np.array([[C_var,C_cov],[C_cov,C_var]])
    data = r.multivariate_normal([a,a],C, size = 100)
    data2 = r.multivariate_normal([b,b],C, size = 100)
    data = np.vstack((data, data2))
    clust = KMeans(n_clusters=2).fit(data)
    plt.scatter(data[:,0], data[:,1], c = clust.predict(data), s =30)
    plt.scatter([a,b],[a,b], s = 80, c = 'green')
    plt.show()
    return


site_mark_dict = {'CB':'o', 'LA':'^', 'MV':'*'}
markers = [site_mark_dict[i] for i in list(data_full.ix[morph_subset.index, 'site'])]
morph_clust_data = np.vstack((m_score[0], m_score[1])).T
clust = KMeans(n_clusters = 2).fit(morph_clust_data)
colors = clust.predict(morph_clust_data)
colors = ['blue' if c == 0 else 'red' for c in colors]
#sizes = env.ix[morph_subset.index, 'alt']
sizes = e_score[0]
sizes =  (sizes - sizes.min())/max(sizes-sizes.min())
sizes = list(80 + 240*sizes)
clust_fig = plt.figure()
clust_fig.suptitle('K-means cluster analysis of morpho data, plotted against morphological PCs\nshape=site (CB:circle, LA:triangle, MV:star)      size = env PC1', size = 20)
ax = clust_fig.add_subplot(111)
for i in range(morph_clust_data.shape[0]):
    ax.scatter(morph_clust_data[i,0], morph_clust_data[i,1], color = colors[i], marker = markers[i], s = sizes[i])
    #plt.text(morph_clust_data[i,0], morph_clust_data[i,1]+0.2, str(env.ix[morph_subset.index[i], 'alt']))

ax.xaxis.set_label_text('morph PC1', size = 15)
ax.yaxis.set_label_text('morph PC2', size = 15)



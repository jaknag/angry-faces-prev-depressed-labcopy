
# Preparations
def add_values_in_dict(dict, key, list_of_values):
    ''' Append multiple values to a key in
        the given dictionary '''
    if key not in dict:
        dict[key] = list_of_values
    else:
        dict[key].extend(list_of_values)
    return dict

def compare_two_proportions(n1,n2,p1,p2,alpha=0.05):
    '''Compare if there is a statistically significant difference between two proportions
        n1 - number of observations in first group
        p1 - number of specific instances in first group
        e.g. 12 students, 5 are female; n1 = 12, p1 = 5
    '''
    p_bar1 = p1 / n1
    p_bar2 = p2 / n2
    p_bar = (p1+p2) / (n1+n2)
    std_error = np.sqrt(p_bar*(1-p_bar)*(1/n1+1/n2))  # standard error
    test_statistic = (p_bar1 - p_bar2)/std_error      # test statistic
    p_value = 2*stats.norm.sf(abs(test_statistic))

    if p_value<=alpha:
        print(f'Reject H0: z-statistic={test_statistic:.3f} p-value={p_value:.4f}')

    if p_value>alpha:
        print(f'Failed to reject H0: z-statistic={test_statistic:.3f} p-value={p_value:.4f}')

    return


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import scipy.stats as stats
import scipy.stats
import statistics
import math
import pingouin as pg

plt.close('all')
beh_path = ('/Users/nagrodzkij/data/angry/')
input_path = ('/Users/nagrodzkij/data/angry/input/')
output_path = ('/Users/nagrodzkij/data/angry/output/')

FaceScores = pd.read_csv(input_path+'/demog/FaceScores.csv',index_col=[0])

prev_depressed = pd.read_csv(input_path+'/demog/previously_depressed.csv', names=['0'])
list_depressed = prev_depressed['0'].to_list()

demog = pd.read_csv(output_path+'table_demog_withedu_byccid.csv')

# Prepare trial data (RT and accuracy)
trialdata = pd.read_csv(output_path+'data_emoFace_excl.csv')
trialdata['correct']=np.where(trialdata['response'] == trialdata['stim_col'], 1, 0)
trialdata_bysubj = trialdata.groupby("subj_idx")['rt'].median().reset_index()
trialdata_bysubj_accuracy =trialdata.groupby('subj_idx')['correct'].mean().reset_index()
trialdata_bysubj = trialdata_bysubj.merge(trialdata_bysubj_accuracy,how='left',on='subj_idx')

# COMPARE THE TWO GROUPS - MANN WHITNEY U VS T-TEST

# Exclude participant who is currently depressed
list_all = list(demog.ccid)
list_nondepressed = list(set(list_all) - set(list_depressed))

demog_depressed = demog[demog['ccid'].isin(list_depressed)]
demog_nondepressed = demog[demog['ccid'].isin(list_nondepressed)]

demog = demog[demog.ccid != 520055]

# RT

trialdata_depressed = trialdata_bysubj[trialdata_bysubj['subj_idx'].isin(list_depressed)]
trialdata_nondepressed = trialdata_bysubj[trialdata_bysubj['subj_idx'].isin(list_nondepressed)]

print('Testing difference in RT - choose the correct test')
data_group1 = trialdata_depressed['rt']
data_group2 = trialdata_nondepressed['rt']

stat1, p1 = stats.shapiro(data_group1)
stat2, p2 = stats.shapiro(data_group2)

print(p1)
print(p2)


if p1 and p2 > 0.05:
    print('Normally distributed, therefore performing t-test')
    tstat = pg.ttest(data_group1,data_group2,correction=True)
    print('Variance: ',np.var(data_group1), np.var(data_group2))
    print('Mean depressed: ' + str(statistics.mean(data_group1)) + 'Mean non-depressed:' + str(statistics.mean(data_group2)))
    print('Std depressed: '+ str(math.sqrt(np.var(data_group1))) + 'Std non-depressed: ' + str(math.sqrt(np.var(data_group2))))
    print(tstat)
    print('')

else:
    print('Not normally distributed, performing Mann Whitney U test')
    print('Median depressed: ' + str(statistics.median(data_group1)))
    q3_depressed, q1_depressed = np.percentile(data_group1, [75 ,25])
    print('IQR depressed: ' + str(q1_depressed) + ', ' + str(q3_depressed))

    print(' ')

    print('Median non-depressed:' + str(statistics.median(data_group2)))
    q3_nondepressed, q1_nondepressed = np.percentile(data_group2, [75 ,25])
    print('IQR non-depressed: ' + str(q1_nondepressed) + ', ' + str(q3_nondepressed))
    res = scipy.stats.mannwhitneyu(data_group1, data_group2)
    print(res)
    print('')

    data_total = pd.concat([data_group1,data_group2],axis=0)

    print('Median all: ' + str(statistics.median(data_total)))
    q3_total,q1_total = np.percentile(data_total, [75 ,25])
    print('IQR non-depressed: ' + str(q1_total) + ', ' + str(q3_total))
    print('')

# ACCURACY

print('Testing difference in accuracy - choose the correct test')
data_group1 = trialdata_depressed['correct']
data_group2 = trialdata_nondepressed['correct']

stat1, p1 = stats.shapiro(data_group1)
stat2, p2 = stats.shapiro(data_group2)

print(p1)
print(p2)


if p1 and p2 > 0.05:
    print('Normally distributed, therefore performing t-test')
    tstat = pg.ttest(data_group1,data_group2,correction=True)
    print('Variance: ',np.var(data_group1), np.var(data_group2))
    print('Mean depressed: ' + str(statistics.mean(data_group1)) + 'Mean non-depressed:' + str(statistics.mean(data_group2)))
    print('Std depressed: '+ str(math.sqrt(np.var(data_group1))) + 'Std non-depressed: ' + str(math.sqrt(np.var(data_group2))))
    print(tstat)
    print('')

else:
    print('Not normally distributed, performing Mann Whitney U test')
    print('Median depressed: ' + str(statistics.median(data_group1)))
    q3_depressed, q1_depressed = np.percentile(data_group1, [75 ,25])
    print('IQR depressed: ' + str(q1_depressed) + ', ' + str(q3_depressed))

    print(' ')

    print('Median non-depressed:' + str(statistics.median(data_group2)))
    q3_nondepressed, q1_nondepressed = np.percentile(data_group2, [75 ,25])
    print('IQR non-depressed: ' + str(q1_nondepressed) + ', ' + str(q3_nondepressed))
    res = scipy.stats.mannwhitneyu(data_group1, data_group2)
    print(res)
    print('')

    data_total = pd.concat([data_group1,data_group2],axis=0)

    print('Median all: ' + str(statistics.median(data_total)))
    q3_total,q1_total = np.percentile(data_total, [75 ,25])
    print('IQR non-depressed: ' + str(q1_total) + ', ' + str(q3_total))
    print('')


# SEX
number_depressed = len(demog_depressed)
number_nondepressed = len(demog_nondepressed)

female_depressed = sum(demog_depressed['sex']=='F')
female_nondepressed = sum(demog_nondepressed['sex']=='F')

male_depressed = sum(demog_depressed['sex']=='M')
male_nondepressed = sum(demog_nondepressed['sex']=='M')

print('Testing if difference in proportion of females - chi squared test')
data_for_test = [[female_depressed,male_depressed],[female_nondepressed,male_nondepressed]]
stat,p,dof,expected = stats.chi2_contingency(data_for_test)
alpha = 0.05
print('X2 statistic is ' +str(stat))
print("p value is " + str(p))
print("dof is "+str(dof))
N = female_depressed + female_nondepressed + male_depressed + male_nondepressed
num_female = female_depressed + female_nondepressed
prop_female = num_female / N
print("N is "+str(N))
print("number female is "+str(num_female))
print("proportion of female is " +str(prop_female))
print("Number female depressed is " + str(female_depressed) + " which is " +str(female_depressed/(female_depressed+male_depressed)) +" of all depressed")
print("Number female non-depressed is " + str(female_nondepressed) + " which is " +str(female_nondepressed/(female_nondepressed+male_nondepressed)) +" of all nondepressed")
print('')
if p <= alpha:
    print('Dependent (reject H0)')
else:
    print('Independent (H0 holds true)')

print('')

#########################
# HANDEDNESS
rh_depressed = sum(demog_depressed['handedness']>0)
rh_nondepressed = sum(demog_nondepressed['handedness']>0)

lh_depressed = sum(demog_depressed['handedness']<=0)
lh_nondepressed = sum(demog_nondepressed['handedness']<=0)


print('Testing if difference in proportion of right handed - chi squared test')
data_for_test = [[rh_depressed,lh_depressed],[rh_nondepressed,lh_nondepressed]]
stat,p,dof,expected = stats.chi2_contingency(data_for_test)
alpha = 0.05
print('X2 statistic is ' +str(stat))
print("p value is " + str(p))
print("dof is "+str(dof))
if p <= alpha:
    print('Dependent (reject H0)')
else:
    print('Independent (H0 holds true)')

print('')

#########################
# EDUCATION
dict_for_test = {}
edu_list = ['None','GCSE','A-Level','Uni']
labels = [list_depressed, list_nondepressed, 'Total']
groups = ['Dep','Not dep']
for i in range(0,2):
    filtered = demog[demog['ccid'].isin(labels[i])]
    number = len(filtered)

    for k in range(0,4):
        edu = sum(filtered['highest_qualification']==k)
        edu_prop = edu/number
        demog_table_dict = add_values_in_dict(dict_for_test, 'Edu', [edu_list[k]])
        demog_table_dict = add_values_in_dict(dict_for_test, 'Group', [groups[i]])
        demog_table_dict = add_values_in_dict(dict_for_test, 'Number', [edu])

df_for_test = pd.DataFrame(data=demog_table_dict)

print('Testing if level of education dependent on group - chi squared test')
data_for_test = [[2,1,3,13],[4,7,10,94]]
stat,p,dof,expected = stats.chi2_contingency(data_for_test)
alpha = 0.05
print('X2 statistic is ' +str(stat))
print("p value is " + str(p))
print("dof is "+str(dof))
if p <= alpha:
    print('Dependent (reject H0)')
else:
    print('Independent (H0 holds true)')

print('')

####################
# ACER
print('Testing difference in ACER - choose the correct test')
data_group1 = demog_depressed['acer']
data_group2 = demog_nondepressed['acer']

stat1, p1 = stats.shapiro(data_group1)
stat2, p2 = stats.shapiro(data_group2)

print(p1)
print(p2)


if p1 and p2 > 0.05:
    print('Normally distributed, therefore performing t-test')
    tstat = pg.ttest(data_group1,data_group2,correction=True)
    print('Variance: ',np.var(data_group1), np.var(data_group2))
    print('Mean depressed: ' + str(statistics.mean(data_group1)) + 'Mean non-depressed:' + str(statistics.mean(data_group2)))
    print('Std depressed: '+ str(math.sqrt(np.var(data_group1))) + 'Std non-depressed: ' + str(math.sqrt(np.var(data_group2))))
    print(tstat)
    print('')

else:
    print('Not normally distributed, performing Mann Whitney U test')
    print('Median depressed: ' + str(statistics.median(data_group1)))
    q3_depressed, q1_depressed = np.percentile(data_group1, [75 ,25])
    print('IQR depressed: ' + str(q1_depressed) + ', ' + str(q3_depressed))

    print(' ')

    print('Median non-depressed:' + str(statistics.median(data_group2)))
    q3_nondepressed, q1_nondepressed = np.percentile(data_group2, [75 ,25])
    print('IQR non-depressed: ' + str(q1_nondepressed) + ', ' + str(q3_nondepressed))
    res = scipy.stats.mannwhitneyu(data_group1, data_group2)
    print(res)
    print('')

    data_total = pd.concat([data_group1,data_group2],axis=0)

    print('Median all: ' + str(statistics.median(data_total)))
    q3_total,q1_total = np.percentile(data_total, [75 ,25])
    print('IQR non-depressed: ' + str(q1_total) + ', ' + str(q3_total))
    print('')


###################
# AGE
print('Testing difference in age - choose the correct test')
data_group1 = demog_depressed['age']
data_group2 = demog_nondepressed['age']
stat1, p1 = stats.shapiro(data_group1)
stat2, p2 = stats.shapiro(data_group2)

print(p1)
print(p2)

if p1 and p2 > 0.05:
    print('Normally distributed, choose the correct test')
    tstat = pg.ttest(data_group1,data_group2,correction=True)
    print('Variance: ',np.var(data_group1), np.var(data_group2))
    print('Mean depressed: ' + str(statistics.mean(data_group1)) + 'Mean non-depressed:' + str(statistics.mean(data_group2)))
    print('Std depressed: '+ str(math.sqrt(np.var(data_group1))) + 'Std non-depressed: ' + str(math.sqrt(np.var(data_group2))))
    print(tstat)
    print('')

else:
    print('Not normally distributed, performing Mann Whitney U test')
    print('Median depressed: ' + str(statistics.median(data_group1)))
    q3_depressed, q1_depressed = np.percentile(data_group1, [75 ,25])
    print('IQR depressed: ' + str(q1_depressed) + ', ' + str(q3_depressed))

    print(' ')

    print('Median non-depressed:' + str(statistics.median(data_group2)))
    q3_nondepressed, q1_nondepressed = np.percentile(data_group2, [75 ,25])
    print('IQR non-depressed: ' + str(q1_nondepressed) + ', ' + str(q3_nondepressed))
    res = scipy.stats.mannwhitneyu(data_group1, data_group2)
    print(res)
    print('')
    data_total = pd.concat([data_group1,data_group2],axis=0)

    print('Median all: ' + str(statistics.median(data_total)))
    q3_total,q1_total = np.percentile(data_total, [75 ,25])
    print('IQR non-depressed: ' + str(q1_total) + ', ' + str(q3_total))
    print('')

###################
# HADS depression
print('Testing difference in HADS depression - choose the correct test')
data_group1 = demog_depressed['hads_depression']
data_group2 = demog_nondepressed['hads_depression']
stat1, p1 = stats.shapiro(data_group1)
stat2, p2 = stats.shapiro(data_group2)

print(p1)
print(p2)

if p1 and p2 > 0.05:
    print('Normally distributed, therefore performing t-test')
    tstat = pg.ttest(data_group1,data_group2,correction=True)
    print('Variance: ',np.var(data_group1), np.var(data_group2))
    print('Mean depressed: ' + str(statistics.mean(data_group1)) + 'Mean non-depressed:' + str(statistics.mean(data_group2)))
    print('Std depressed: '+ str(math.sqrt(np.var(data_group1))) + 'Std non-depressed: ' + str(math.sqrt(np.var(data_group2))))
    print(tstat)
    print('')

else:
    print('Not normally distributed, performing Mann Whitney U test')
    print('Median depressed: ' + str(statistics.median(data_group1)))
    q3_depressed, q1_depressed = np.percentile(data_group1, [75 ,25])
    print('IQR depressed: ' + str(q1_depressed) + ', ' + str(q3_depressed))

    print(' ')

    print('Median non-depressed:' + str(statistics.median(data_group2)))
    q3_nondepressed, q1_nondepressed = np.percentile(data_group2, [75 ,25])
    print('IQR non-depressed: ' + str(q1_nondepressed) + ', ' + str(q3_nondepressed))
    res = scipy.stats.mannwhitneyu(data_group1, data_group2)
    print(res)
    print('')
    data_total = pd.concat([data_group1,data_group2],axis=0)

    print('Median all: ' + str(statistics.median(data_total)))
    q3_total,q1_total = np.percentile(data_total, [75 ,25])
    print('IQR non-depressed: ' + str(q1_total) + ', ' + str(q3_total))
    print('')

###################
# Benton
print('Testing difference in Benton scores - choose the correct test')
FaceScores_depressed = FaceScores.loc[list_depressed,:]
data_group1 = FaceScores_depressed['benton']
FaceScores_nondepressed = FaceScores.loc[list_nondepressed,:]
data_group2 = FaceScores_nondepressed['benton']
stat1, p1 = stats.shapiro(data_group1)
stat2, p2 = stats.shapiro(data_group2)


print(p1)
print(p2)


if p1 and p2 > 0.05:
    print('Normally distributed, therefore performing t-test')
    tstat = pg.ttest(data_group1,data_group2,correction=True)
    print('Variance: ',np.var(data_group1), np.var(data_group2))
    print('Mean depressed: ' + str(statistics.mean(data_group1)) + 'Mean non-depressed:' + str(statistics.mean(data_group2)))
    print('Std depressed: '+ str(math.sqrt(np.var(data_group1))) + 'Std non-depressed: ' + str(math.sqrt(np.var(data_group2))))
    print(tstat)
    print('')

else:
    print('Not normally distributed, performing Mann Whitney U test')
    print('Median depressed: ' + str(statistics.median(data_group1)))
    q3_depressed, q1_depressed = np.percentile(data_group1, [75 ,25])
    print('IQR depressed: ' + str(q1_depressed) + ', ' + str(q3_depressed))

    print(' ')

    print('Median non-depressed:' + str(statistics.median(data_group2)))
    q3_nondepressed, q1_nondepressed = np.percentile(data_group2, [75 ,25])
    print('IQR non-depressed: ' + str(q1_nondepressed) + ', ' + str(q3_nondepressed))
    res = scipy.stats.mannwhitneyu(data_group1, data_group2)
    print(res)
    print('')
    data_total = pd.concat([data_group1,data_group2],axis=0)

    print('Median all: ' + str(statistics.median(data_total)))
    q3_total,q1_total = np.percentile(data_total, [75 ,25])
    print('IQR non-depressed: ' + str(q1_total) + ', ' + str(q3_total))
    print('')